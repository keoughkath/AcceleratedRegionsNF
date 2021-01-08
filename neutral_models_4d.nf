VERSION = "0.0.0"

// parse inputs
// use species-filtered mafs if you want

initMafChannel = Channel.fromPath( "${params.maf_path}chr*.maf.gz" ).filter( ~/.*chr.{1,2}\.maf\.gz/ )
speciesFile = file(params.species)
speciesList = speciesFile.readLines().join(",")
initTree = file(params.init_tree)

/*
* extract the species of interest from the overall MAF
*/

process extractSpecies {
	tag "Extracting species of interest from the MAF for ${initMaf}"

	// publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(species_list_file) from speciesFile
	file(initMaf) from initMafChannel

	output:
	file("*_species.maf") into speciesMaf

	script:
	//
	// UCSC MafSpeciesSubset
	//
	chrom = initMaf.simpleName
	"""
	mafSpeciesSubset ${params.maf_path}${initMaf} ${species_list_file} ${chrom}_species.maf
	"""

}

/*
* extract the 4d codons from the species-extracted MAF
*/

process extract4dCodons {
	tag "Extracting 4D codons for ${chrom}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(maf) from speciesMaf

	output:
	file("*_4dcodons.ss") into fourDCodons

	script:
	// 
	// msa_view from PHAST
	//
	chrom = maf.simpleName.split('_')[0]
	"""
	${params.phast_path}./msa_view --in-format MAF --out-format SS --4d --features ${params.gtf_path}/${chrom}_refseq_genes_hg38.gtf ${maf} > ${chrom}_4dcodons.ss
	"""

}


/*
* extract 4d sites from 4d codons
*/

process extract4dSites{
	tag "Extracting 4D sites for ${chrom}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(codons) from fourDCodons

	output:
	file("*_4dsites.ss") into fourDSites

	script:
	//
	// msa_view from PHAST
	//
	chrom = codons.simpleName.split('_')[0]
	"""
	${params.phast_path}./msa_view --in-format SS --out-format SS --tuple-size 1 ${codons} > ${chrom}_4dsites.ss
	"""

}


// copy fourDSites into 2 channels and filter by autosomes and non-autosomes

fourDSites.into { fourDSitesOne; fourDSitesTwo}

fourDSitesAuto = fourDSitesOne.filter( ~/.*chr\d.*/ )
fourDSitesNonAuto = fourDSitesTwo.filter( ~/.*chr[X Y].*/ )

/*
* aggregate 4d sites from autosomes into a single file
*/

process aggregateAutosome4dSites {
	tag "Aggregating 4d sites from autosomes"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	val(species_list) from speciesList
	val(auto_ss_files) from fourDSitesAuto.collect()

	output:
	file("autosomes_4dsites.ss") into autoFourDSitesAgg

	script:
	//
	// msa_view from PHAST
	//
	ss_file_list = auto_ss_files.collect{"$it"}.join(",")
	"""
	${params.phast_path}./msa_view --unordered-ss --out-format SS --aggregate ${species_list} ${ss_file_list} > autosomes_4dsites.ss
	"""

}

/*
* prune the species tree to reflect only the species of interest
*/

process pruneTree {
	tag "Pruning the species tree"

	// publishDir params.outdir, mode: "copy", overwrite: false

	input:
	file(tree) from initTree
	val(species_list) from speciesList

	output:
	file("pruned_tree.nh") into prunedTree

	script:
	// 
	// tree_doctor from PHAST
	//
	"""
	${params.phast_path}./tree_doctor -P ${species_list} ${tree} > pruned_tree.nh
	"""
}


/*
* generate neutral model using 4d sites
* autosomes model separate from sex chromosomes
*/

process neutralModelAutosomes {
	tag "Generating neutral model for autosomes"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input: 
	file(tree) from prunedTree
	file(four_d_sites) from autoFourDSitesAgg

	output:
	file('autosomes_neutral_model_prelim.mod') into autoNeutralModelPrelim

	script:
	//
	// phyloFit from PHAST
	//
	"""
	${params.phast_path}./phyloFit --tree ${tree} --out-root autosomes_neutral_model_prelim --msa-format SS --EM --seed ${params.random_seed} --log autosomes.log ${four_d_sites}
	"""
}



process neutralModelNonAutosomes {
	tag "Generating neutral models for non-autosomes"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(tree) from prunedTree
	file(four_d_sites) from fourDSitesNonAuto

	output:
	file('*_neutral_model_prelim.mod') into nonAutoNeutralModelPrelim

	script:
	//
	// phyloFit from PHAST
	//
	chrom = four_d_sites.simpleName.split('_')[0]
	"""
	${params.phast_path}./phyloFit --tree ${tree} --out-root ${chrom}_neutral_model_prelim --msa-format SS --EM --seed ${params.random_seed} --log ${chrom}.log ${four_d_sites}
	"""

}


/*
* change base frequencies in neutral models to reflect the whole genome
* rather than just the 4D sites
*/

process modifyAutoBaseFrequencies {
	tag "Modifying base frequencies of neutral models for autosomes"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(prelim_neutral_model) from autoNeutralModelPrelim

	output:
	file("autosomes.mod") into autoNeutralModel

	script:
	//
	// modFreqs from PHAST
	//
	"""
	${params.phast_path}./modFreqs ${prelim_neutral_model} 0.3 0.2 0.2 0.3 > autosomes.mod
	"""
}

autoNeutralModel.into { autoNeutralModelOne; autoNeutralModelTwo }


process modifyNonAutoBaseFrequencies {
	tag "Modifying base frequencies of neutral models for non-autosome ${chrom} "

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(prelim_neutral_model) from nonAutoNeutralModelPrelim

	output:
	file("*.mod") into nonAutoNeutralModels

	script:
	//
	// modFreqs from PHAST
	//
	chrom = prelim_neutral_model.simpleName.split('_')[0]
	"""
	${params.phast_path}./modFreqs ${prelim_neutral_model} 0.2 0.3 0.3 0.2 > ${chrom}.mod
	"""
}