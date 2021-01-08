VERSION = "0.0.2"

// parse inputs

speciesFile = file(params.species)
initTree = file(params.init_tree)
speciesList = speciesFile.readLines().join(",")
chromBedPath = file(params.chrom_bed_path)
syntenyFiles = Channel.fromPath( "${params.synteny_filter_path}*.bed" )
arFilterFile = file(params.ar_filters_path)
// autoNeutralModel = Channel.fromPath( params.auto_neutral_model )
// nonAutoNeutralModels = Channel.fromPath( "${params.nonauto_neutral_model}*.mod" )

// autoNeutralModel.into { autoNeutralModelOne; autoNeutralModelTwo }

/*
* set up channels for each chromosome for parallelization, filter out non-standard chromosome MAFs
*/

initMafChannel = Channel.fromPath( "${params.maf_path}chr*.maf.gz" )//.filter( ~/.*chr.{1,2}\.maf\.gz/ )

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
	fname = initMaf.baseName
	"""
	mafSpeciesSubset ${params.maf_path}${initMaf} ${species_list_file} ${fname}_species.maf
	"""

}

// copy speciesMaf into 5 channels to use in 5 processes

speciesMaf.into { speciesMafOne; speciesMafTwo; speciesMafThree; speciesMafFour; speciesMafFive }

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

process maskSpeciesOfInterest {
	tag "Masking ${params.species_of_interest} for ${fname}."

	// publishDir params.outdir, mode: "copy", overwrite: true

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
	fname = species_maf.baseName
	chrom = species_maf.simpleName
	"""
	${params.phast_path}./maf_parse --features ${chrom_bed_path}/hg38_${chrom}.bed --mask-features ${params.species_of_interest} ${species_maf} > ${fname}_${params.species_of_interest}_masked.maf
	"""
}

// copy maskedMaf into 2 channels and filter for autosomal and non-autosomal

maskedMaf.into { maskedMafOne; maskedMafTwo }

maskedMafsAuto = maskedMafOne.filter( ~/.*chr\d.*/ )
maskedMafsNonAuto = maskedMafTwo.filter( ~/.*chr[X Y M].*/ )

// join non-autosomal neutral model and masked_maf channels for non-autosomes so they are coordinated

// maskedMafsNonAutoChr = maskedMafsNonAuto.map {
// 	file -> tuple(file.simpleName, file)
// }

// nonAutoNeutralModelsChr = nonAutoNeutralModels.filter( ~/.*chr[X Y M].*/ ).map {
// 	file -> tuple(file.simpleName, file)
// }

// nonAutoNeutralModelsChr.into { nonAutoNeutralModelsChrOne; nonAutoNeutralModelsChrTwo }

// nonAutoMaskedNeutral = maskedMafsNonAutoChr.join(nonAutoNeutralModelsChrOne)

// maskedMafsAuto.subscribe { println it }

/*
* run phastCons to identify elements conserved 
*/

process callAutosomalConservedElements {
	tag "Identifying conserved loci on ${fname}."

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(masked_maf) from maskedMafsAuto
	// file(neutral_model) from autoNeutralModelOne

	output:
	file("*_phastCons_unfilt.bed") into autoPhastConsElements

	script:
	//
	// phastCons from PHAST - make sure you use the github version :/
	//
	chrom = masked_maf.simpleName
	fname = masked_maf.baseName
	"""
	${params.phast_path}./phastCons ${masked_maf} ${params.auto_neutral_model} --rho 0.3 --target-coverage ${params.target_coverage} --expected-length ${params.expected_length} --no-post-probs --score --viterbi ${fname}_phastCons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

process callNonAutosomalConservedElements {
	tag "Identifying conserved loci on non-autosomal ${fname}."

	// publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	set file(masked_maf) from maskedMafsNonAuto

	output:
	file("${fname}_phastCons_unfilt.bed") into nonAutoPhastConsElements

	script:
	//
	// phastCons from PHAST - make sure you use the github version :/
	//
	fname = masked_maf.baseName
	chrom = masked_maf.simpleName
	"""
	${params.phast_path}./phastCons ${masked_maf} ${params.nonauto_neutral_model}${chrom}.neutral.mod --rho 0.3 --target-coverage ${params.target_coverage} --expected-length ${params.expected_length} --no-post-probs --score --viterbi ${fname}_phastCons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

// merge autosomal and non-autosomal phastCons channels and filter instances where there aren't any conserved elements

phastCons = autoPhastConsElements.concat(nonAutoPhastConsElements).filter{ it.size()>0 }

// would be good here to merge phastCons from different windows on same chromosome...

/*
* filter phastCons by log-odds score
*/

process filterByScore {
	tag "Filtering phastcons from ${fname} by log odds score"

	input:
	file(phastcons) from phastCons

	output:
	file("*_score_filtered.bed") into phastConsScoreFiltered

	script:
	chrom = phastcons.simpleName
	fname = phastcons.baseName
	"""
	python ${baseDir}/bin/score_filter_phastcons.py ${phastcons} ${fname}_phastcons_score_filtered.bed
	"""
}

/*
* grab phastCons form regions that are syntenic
* with user-identified species
*/

process syntenyFilterPhastCons {
	tag "Synteny filtering phastCons from ${fname}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsScoreFiltered.filter{ it.size()>0 }
	file(synteny_file_list) from syntenyFiles.collect()

	output:
	file("*_phastcons_synteny_filtered.bed") into phastConsSyntenyFiltered

	script:
	chrom = phastcons.simpleName
	fname = phastcons.baseName
	"""
	bedops -i ${phastcons} ${synteny_file_list} > ${fname}_phastcons_synteny_filtered.bed
	"""
}

/*
* filter phastCons for elements that could confound the acceleration analysis
* e.g. duplicated regions, pseudogenes, etc.
*/

process phastconsDupFilter {
	tag "Filtering phastCons from ${fname} for repeats, duplicated regions, pseudogenes and selfChain"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsSyntenyFiltered.filter{ it.size()>0 }
	file(ar_filters) from arFilterFile

	output:
	file("*_phastcons_dups_filtered.bed") into phastConsFiltered

	script:
	chrom = phastcons.simpleName
	fname = phastcons.baseName
	"""
	bedtools subtract -a ${phastcons} -b ${ar_filters} | bedtools sort -i - | bedtools merge -d 0 -i - > ${fname}_phastcons_dups_filtered.bed
	"""
}


/*
* filter phastCons elements for minimum size of element
*/

process phastconsSizeFilter {
	tag "Filtering phastCons from ${fname} for size > 50bp"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsFiltered.filter{ it.size()>0 }

	output:
	file("*_phastcons_filtered.bed") into phastConsSizeFilteredPreID

	script:
	chrom = phastcons
	fname = phastcons.baseName
	"""
	awk '(\$3-\$2) >= 50' ${phastcons} > ${fname}_phastcons_filtered.bed
	"""
}


//* add ID to phastCons, necessary for extracting MAFs


process idToPhastcons {
	tag "Adding ID to phastCons from ${fname}"

	// publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastcons) from phastConsSizeFilteredPreID.filter{ it.size()>0 }

	output:
	file("*_phastcons_filtered_id.bed") into phastConsSizeFiltered

	script:
	chrom = phastcons.simpleName
	fname = phastcons.baseName.split(/.maf/)[0]
	"""
	awk -F'\t' -v OFS='\t' '{ \$(NF+1)=\$1"_"NR} 1' ${phastcons} > ${fname}_phastcons_filtered_id.bed
	"""

}

// /*
// * split phastCons elements into multiple files on which to run phyloP in order to speed the process up 
// */

// process splitPhastCons {
// 	tag "Splitting up phastCons from ${fname} to speed up phyloP"

// 	input:
// 	file(phastcons) from phastConsSizeFiltered.filter{ it.size()>0 } // filters out any instances where filtering has removed all phastCons, e.g. chrY

// 	output:
// 	file("*.bed") into phastConsSplit

// 	script:
// 	chrom = phastcons.simpleName
// 	fname = phastcons.baseName
// 	"""
// 	bedtools split -i ${phastcons} -n ${params.n_phylop_chunks} -p ${fname}_ -a simple
// 	"""
// }

// join neutral model, split phastCons elements and unmasked MAF channels so they are coordinated

phastConsSplit = phastConsSizeFiltered.filter{ it.size()>0 } // filters out any instances where filtering has removed all phastCons, e.g. chrY

phastConsSplit.into { phastConsSplitOne; phastConsSplitTwo }

phastConsSplitAuto = phastConsSplitOne.flatten().filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.baseName.split(/_phastcons/)[0], file)
}

speciesMafChromAuto = speciesMafThree.filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.baseName.split(/.maf/)[0], file)
}

phastConsSplitAutoMaf = phastConsSplitAuto.combine(speciesMafChromAuto, by: 0)//.combine(autoNeutralModelTwo)

// phastConsSplitAutoMaf.subscribe { println it }


// non-autosomes


phastConsSplitNonAuto = phastConsSplitTwo.flatten().filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.baseName.split(/_phastcons/)[0], file)
}

speciesMafChromNonAuto = speciesMafFour.filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.baseName.split(/.maf/)[0], file)
}

phastConsSplitNonAutoMaf = phastConsSplitNonAuto.combine(speciesMafChromNonAuto, by: 0)//.combine(nonAutoNeutralModelsChrTwo, by: 0)

// phastConsSplitNonAutoMaf.subscribe { println it }


/*
* run phyloP to identify accelerated elements, autosomes
*/

process accRegionsAutosomal {
	tag "Identifying accelerated regions on ${window}"

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(window), file(phastcons), file(species_maf) from phastConsSplitAutoMaf

	output:
	stdout into phyloPResultsAuto

	script:
	//
	// phyloP from PHAST
	//
	chrom = phastcons.simpleName
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -g ${params.auto_neutral_model} ${species_maf}
	"""
}

/*
* run phyloP to identify accelerated elements, non-autosomes
*/


process accRegionsNonAutosomal {
	tag "Identifying accelerated regions on ${window}"

	errorStrategy 'retry'
	maxRetries 3

	input:
	set val(window), file(phastcons), file(species_maf) from phastConsSplitNonAutoMaf

	output:
	stdout into phyloPResultsNonAuto

	script:
	//
	// phyloP from PHAST
	//
	chrom = phastcons.simpleName
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -g ${params.nonauto_neutral_model}${chrom}.neutral.mod ${species_maf}
	"""

}

// concatenate accelerated regions

phastConsScored = phyloPResultsAuto.concat(phyloPResultsNonAuto).collectFile(name: 'phyloP_results.txt')

/*
* multiple test correction with Benjamini-Hochberg methodology
*/


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
	python ${baseDir}/bin/multiple_test_correction_filtering.py ${scored_phastcons} ${params.max_p} final_ARs_${params.random_seed}.bed scored_phastCons_${params.random_seed}.txt
	"""

}



