VERSION = "0.0.0"

// parse inputs

speciesFile = file(params.species)
initTree = file(params.init_tree)
speciesList = speciesFile.readLines().join(",")
chromBedPath = file(params.chrom_bed_path)
syntenyFilterFileDog = file(params.synteny_filter_dog)
syntenyFilterFileRhesus = file(params.synteny_filter_rhesus)
syntenyFilterFileMouse = file(params.synteny_filter_mouse)
arFilterFile = file(params.ar_filters_path)

/*
* set up channels for each chromosome for parallelization
*/

initMafChannel = Channel.fromPath( "${params.maf_path}chr*.maf.gz" )

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

// copy speciesMaf into 4 channels to use in 4 processes

speciesMaf.into { speciesMafOne; speciesMafTwo; speciesMafThree; speciesMafFour }

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
* extract the 4d codons from the species-extracted MAF
*/

process extract4dCodons {
	tag "Extracting 4D codons for ${chrom}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(maf) from speciesMafOne

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
fourDSitesNonAuto = fourDSitesTwo.filter( ~/.*chr[X Y M].*/ )

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

	// publishDir params.outdir, mode: "copy", overwrite: true

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

	// publishDir params.outdir, mode: "copy", overwrite: true

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


/*
* mask species of interest (e.g. human) 
*/ 

process maskSpeciesOfInterest {
	tag "Masking ${params.species_of_interest} for ${chrom}."

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(species_maf) from speciesMafTwo
	val(chrom_bed_path) from chromBedPath

	output:
	file("*_${params.species_of_interest}_masked.maf") into maskedMaf

	script:
	//
	// maf_parse from PHAST
	//
	chrom = species_maf.simpleName.split('_')[0]
	"""
	${params.phast_path}./maf_parse --features ${chrom_bed_path}/hg38_${chrom}.bed --mask-features ${params.species_of_interest} ${species_maf} > ${chrom}_${params.species_of_interest}_masked.maf
	"""
}

// copy maskedMaf into 2 channels and filter for autosomal and non-autosomal

maskedMaf.into { maskedMafOne; maskedMafTwo }

maskedMafsAuto = maskedMafOne.filter( ~/.*chr\d.*/ )
maskedMafsNonAuto = maskedMafTwo.filter( ~/.*chr[X Y M].*/ )

// join non-autosomal neutral model and masked_maf channels for non-autosomes so they are coordinated

maskedMafsNonAutoChr = maskedMafsNonAuto.map {
	file -> tuple(file.simpleName.split('_')[0], file)
}

nonAutoNeutralModelsChr = nonAutoNeutralModels.map {
	file -> tuple(file.simpleName, file)
}

nonAutoNeutralModelsChr.into { nonAutoNeutralModelsChrOne; nonAutoNeutralModelsChrTwo }

nonAutoMaskedNeutral = maskedMafsNonAutoChr.join(nonAutoNeutralModelsChrOne)


/*
* run phastCons to identify elements conserved 
*/

process callAutosomalConservedElements {
	tag "Identifying conserved loci on ${chrom}."

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(masked_maf) from maskedMafsAuto
	file(neutral_model) from autoNeutralModelOne

	output:
	file("*_phastCons_unfilt.bed") into autoPhastConsElements

	script:
	//
	// phastCons from PHAST - make sure you use the github version :/
	//
	chrom = masked_maf.simpleName.split('_')[0]
	"""
	${params.phast_path}./phastCons ${masked_maf} ${neutral_model} --rho 0.3 --target-coverage 0.3 --expected-length 45 --no-post-probs --viterbi ${chrom}_phastCons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

process callNonAutosomalConservedElements {
	tag "Identifying conserved loci on non-autosomal ${chrom}."

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(chrom), file(masked_maf), file(neutral_model) from nonAutoMaskedNeutral

	output:
	file("${chrom}_phastCons_unfilt.bed") into nonAutoPhastConsElements

	script:
	//
	// phastCons from PHAST - make sure you use the github version :/
	//
	"""
	${params.phast_path}./phastCons ${masked_maf} ${neutral_model} --rho 0.3 --target-coverage 0.3 --expected-length 45 --no-post-probs --viterbi ${chrom}_phastCons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

// merge autosomal and non-autosomal phastCons channels

phastCons = autoPhastConsElements.concat(nonAutoPhastConsElements)

/*
* grab phastCons form regions that are syntenic
* with user-identified species
*/

process syntenyFilterPhastCons {
	tag "Synteny filtering phastCons from ${chrom}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastCons
	file(synteny_filter_dog) from syntenyFilterFileDog
	file(synteny_filter_rhesus) from syntenyFilterFileRhesus
	file(synteny_filter_mouse) from syntenyFilterFileMouse

	output:
	file("*_phastcons_synteny_filtered.bed") into phastConsSyntenyFiltered

	script:
	chrom = phastcons.simpleName.split('_')[0]
	"""
	bedtools intersect -a ${phastcons} -b ${synteny_filter_dog} -sorted | bedtools intersect -a - -b ${synteny_filter_rhesus} | bedtools intersect -a - -b ${synteny_filter_mouse} | bedtools sort -i - > ${chrom}_phastcons_synteny_filtered.bed
	"""
}

/*
* filter phastCons for elements that could confound the acceleration analysis
* e.g. duplicated regions, pseudogenes, etc.
*/

process phastconsDupFilter {
	tag "Filtering phastCons from ${chrom} for repeats, duplicated regions, pseudogenes and selfChain"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsSyntenyFiltered
	file(ar_filters) from arFilterFile

	output:
	file("*_phastcons_dups_filtered.bed") into phastConsFiltered

	script:
	chrom = phastcons.simpleName.split('_')[0]
	"""
	bedtools subtract -a ${phastcons} -b ${ar_filters} | bedtools sort -i - | bedtools merge -d 0 -i - > ${chrom}_phastcons_dups_filtered.bed
	"""
}


/*
* filter phastCons elements for minimum size of element
*/

process phastconsSizeFilter {
	tag "Filtering phastCons from ${chrom} for size > 50bp"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsFiltered

	output:
	file("*_phastcons_filtered.bed") into phastConsSizeFiltered

	script:
	chrom = phastcons.simpleName.split('_')[0]
	"""
	awk '(\$3-\$2) >= 50' ${phastcons} > ${chrom}_phastcons_filtered.bed
	"""
}

// phastConsSizeFiltered.subscribe { println it.size() }

/*
* split phastCons elements into multiple files on which to run phyloP in order to speed the process up 
*/

process splitPhastCons {
	tag "Splitting up phastCons from ${chrom} to speed up phyloP"

	input:
	file(phastcons) from phastConsSizeFiltered.filter{ it.size()>0 } // filters out any instances where filtering has removed all phastCons, e.g. chrY

	output:
	file("*.bed") into phastConsSplit

	script:
	chrom = phastcons.simpleName.split('_')[0]
	"""
	bedtools split -i ${phastcons} -n ${params.n_phylop_chunks} -p ${chrom}_ -a simple
	"""
}

// join neutral model, split phastCons elements and masked_maf channels so they are coordinated

phastConsSplit.into { phastConsSplitOne; phastConsSplitTwo }

phastConsSplitAuto = phastConsSplitOne.flatten().filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.simpleName.split('_')[0], file)
}

speciesMafChromAuto = speciesMafThree.filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.simpleName.split('_')[0], file)
}

phastConsSplitAutoMafNeutral = phastConsSplitAuto.combine(speciesMafChromAuto, by: 0).spread(autoNeutralModelTwo)


// non-autosomes


phastConsSplitNonAuto = phastConsSplitTwo.flatten().filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.simpleName.split('_')[0], file)
}

speciesMafChromNonAuto = speciesMafFour.filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.simpleName.split('_')[0], file)
}

phastConsSplitNonAutoMafNeutral = phastConsSplitNonAuto.combine(speciesMafChromNonAuto, by: 0).combine(nonAutoNeutralModelsChrTwo, by: 0)


/*
* run phyloP to identify accelerated elements
*/

process accRegionsAutosomal {
	tag "Identifying accelerated regions on autosomes"

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(chrom), file(phastcons), file(unmasked_maf), file(neutral_model) from phastConsSplitAutoMafNeutral

	output:
	stdout into phyloPResultsAuto

	script:
	//
	// phyloP from PHAST
	//
	id = phastcons.simpleName
	chrom = phastcons.simpleName.split('_')[0]
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -g ${neutral_model} ${unmasked_maf}
	"""
}


process accRegionsNonAutosomal {
	tag "Identifying accelerated regions on non-autosomes"

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(chrom), file(phastcons), file(unmasked_maf), file(neutral_model) from phastConsSplitNonAutoMafNeutral

	output:
	stdout into phyloPResultsNonAuto

	script:
	//
	// phyloP from PHAST
	//
	id = phastcons.simpleName
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -g ${neutral_model} ${unmasked_maf}
	"""

}

// concatenate accelerated regions

phastConsScored = phyloPResultsAuto.concat(phyloPResultsNonAuto).collectFile(name: 'phyloP_results.txt')


// multiple test correction with Benjamini-Hochberg methodology


process multipleTestCorrectionFiltering {
	tag "Multiple test correction and filtering"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(scored_phastcons) from phastConsScored

	output:
	file("final_ARs_${params.random_seed}.bed")
	file("scored_phastCons_${params.random_seed}.txt")

	script:
	"""
	python ${baseDir}/multiple_test_correction_filtering.py ${scored_phastcons} ${params.max_p} final_ARs_${params.random_seed}.bed scored_phastCons_${params.random_seed}.txt
	"""

}



