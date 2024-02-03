VERSION = "0.0.1"

nextflow.enable.dsl = 2

// sample run command
// nextflow run neutral_models_4d.nf -w '20way_nm_workdir/' -profile local -params-file neutral_models_config_20way_example.yml 

/*
* extract the species of interest from the overall MAF
*/

process extractSpecies {
	tag "Extracting species of interest from the MAF for ${initMaf}"

	publishDir params.outdir, mode: "copy", overwrite: false

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(speciesFile)
	file(initMaf)

	output:
	file("*_species.maf")

	script:
	//
	// UCSC MafSpeciesSubset
	//
	chrom = initMaf.simpleName
	"""
	mafSpeciesSubset ${initMaf} ${speciesFile} ${chrom}_species.maf
	"""
}

/*
* extract the 4d codons from the species-extracted MAF
*/

process extract4dCodons {
	tag "Extracting 4D codons for ${chrom}"

	publishDir params.outdir, mode: "copy", overwrite: true

	// errorStrategy 'retry'
	// maxRetries 3

	input:
	file(speciesMaf)
	file(gtf)

	output:
	file("*_4dcodons.ss")

	script:
	// 
	// msa_view from PHAST
	//
	chrom = speciesMaf.simpleName.split('_')[0]
	"""
	grep "^${chrom}\t" ${gtf} | sed 's/^/${params.reference_species}./' > ${chrom}.gtf
	${params.phast_path}./msa_view --in-format MAF --out-format SS --4d --features ${chrom}.gtf ${speciesMaf} > ${chrom}_4dcodons.ss
	"""
}

/*
* extract 4d sites from 4d codons
*/

process extract4dSites{
	tag "Extracting 4D sites for ${chrom}"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(fourDCodons)

	output:
	file("*_4dsites.ss") 

	script:
	//
	// msa_view from PHAST
	//
	chrom = fourDCodons.simpleName.split('_')[0]
	"""
	${params.phast_path}./msa_view --in-format SS --out-format SS --tuple-size 1 ${fourDCodons} > ${chrom}_4dsites.ss
	"""
}

/*
* aggregate 4d sites from autosomes into a single file
*/

process aggregateAutosome4dSites {
	tag "Aggregating 4d sites from autosomes"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	val(speciesList)
	file(fourDSitesAuto)

	output:
	file("autosomes_4dsites.ss")

	script:
	//
	// msa_view from PHAST
	//
	ss_file_list = fourDSitesAuto.collect().join(",")
	"""
	${params.phast_path}./msa_view --unordered-ss --out-format SS --aggregate ${speciesList} ${ss_file_list} > autosomes_4dsites.ss
	"""
}

/*
* prune the species tree to reflect only the species of interest
*/

process pruneTree {
	tag "Pruning the species tree"

	// publishDir params.outdir, mode: "copy", overwrite: false

	input:
	file(initTree)
	val(speciesList)

	output:
	file("pruned_tree.nh")

	script:
	// 
	// tree_doctor from PHAST
	//
	"""
	${params.phast_path}./tree_doctor -P ${speciesList} ${initTree} > pruned_tree.nh
	"""
}


/*
* generate neutral model using 4d sites
* autosomes model separate from sex chromosomes
*/

process neutralModelAutosomes {
	tag "Generating neutral model for autosomes"

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input: 
	file(prunedTree)
	file(autoFourDSitesAgg)

	output:
	file('autosomes_neutral_model_prelim.mod')

	script:
	//
	// phyloFit from PHAST
	//
	"""
	${params.phast_path}./phyloFit --tree ${prunedTree} --out-root autosomes_neutral_model_prelim --msa-format SS --EM --seed ${params.random_seed} --log autosomes.log ${autoFourDSitesAgg}
	"""
}



process neutralModelNonAutosomes {
	tag "Generating neutral models for non-autosomes"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(prunedTree)
	file(fourDSitesNonAuto)

	output:
	file('*_neutral_model_prelim.mod') 

	script:
	//
	// phyloFit from PHAST
	//
	chrom = fourDSitesNonAuto.simpleName.split('_')[0]
	"""
	${params.phast_path}./phyloFit --tree ${prunedTree} --out-root ${chrom}_neutral_model_prelim --msa-format SS --EM --seed ${params.random_seed} --log ${chrom}.log ${fourDSitesNonAuto}
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
	file(autoNeutralModelPrelim)

	output:
	file("autosomes.mod")

	script:
	//
	// modFreqs from PHAST
	//
	"""
	${params.phast_path}./modFreqs ${autoNeutralModelPrelim} 0.3 0.2 0.2 0.3 > autosomes.mod
	"""
}


process modifyNonAutoBaseFrequencies {
	tag "Modifying base frequencies of neutral models for non-autosome ${chrom} "

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(nonAutoNeutralModelPrelim)

	output:
	file("*.mod")

	script:
	//
	// modFreqs from PHAST
	//
	chrom = nonAutoNeutralModelPrelim.simpleName.split('_')[0]
	"""
	${params.phast_path}./modFreqs ${nonAutoNeutralModelPrelim} 0.2 0.3 0.3 0.2 > ${chrom}.mod
	"""
}

// parse inputs
speciesFile = file(params.species)
speciesList = speciesFile.readLines().join(",")
initTree = file(params.init_tree)
gtf = file(params.gtf)

workflow {

	// set up a channel per chromosome
	initMaf = Channel.fromPath( "${params.maf_path}chr*.maf.gz" ).filter( ~/.*chr.{1,2}\.maf\.gz/ )

	// filter MAF to only include species of interest
	extractSpecies(speciesFile, initMaf).set{ speciesMaf }

	// extract the 4d codons from the species-extracted MAF
	extract4dCodons(speciesMaf, gtf).set{ fourDCodons }

	// extract 4d sites
	extract4dSites(fourDCodons).set{ fourDSites }

	// filter by autosomes and non-autosomes
	fourDSitesAuto = fourDSites.filter( ~/.*chr\d.*/ )
	fourDSitesNonAuto = fourDSites.filter( ~/.*chr[X Y].*/ )

	// aggregate 4d sites from autosomes into a single file
	aggregateAutosome4dSites(speciesList, fourDSitesAuto).set{ autoFourDSitesAgg }

	// prune the tree to only represent species of interest
	pruneTree(initTree, speciesList).set{ prunedTree }

	// calculate neutral model for autosomes
	neutralModelAutosomes(prunedTree, autoFourDSitesAgg).set{ autoNeutralModelPrelim }

	// calculate neutral model for non-autosomes
	neutralModelNonAutosomes(prunedTree, fourDSitesNonAuto).set{ nonAutoNeutralModelPrelim }

	// change base frequencies in neutral models to reflect the whole genome
	// rather than just the 4D sites
	modifyAutoBaseFrequencies(autoNeutralModelPrelim).set{ autoNeutralModel }
	modifyNonAutoBaseFrequencies(nonAutoNeutralModelPrelim).set{ nonAutoNeutralModel }


}






