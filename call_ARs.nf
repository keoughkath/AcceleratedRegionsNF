VERSION = "0.0.4"

nextflow.enable.dsl = 2

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
	file("*.maf") 

	script:
	//
	// UCSC MafSpeciesSubset
	//
	chrom = initMaf.simpleName
	fname = initMaf.baseName
	"""
	mafSpeciesSubset ${initMaf} ${speciesFile} ${fname}
	"""

}

process maskSpeciesOfInterest {
	tag "Masking ${params.species_of_interest} for ${fname}."

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(speciesMaf)
	val(chromBedPath)

	output:
	file("*.maf")

	script:
	//
	// maf_parse from PHAST
	//
	fname = speciesMaf.baseName
	chrom = speciesMaf.simpleName
	"""
	${params.phast_path}./maf_parse --features ${chromBedPath}/${chrom}.bed --mask-features ${params.species_of_interest} ${speciesMaf} > ${fname}_masked.maf
	"""
}

// /*
// * run phastCons to identify elements conserved 
// */

process callAutosomalConservedElements {
	tag "Identifying conserved loci on ${fname}."

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(maskedMaf)
	file(autoNeutralModel)

	output:
	file("*_unfilt.bed") 

	script:
	//
	// phastCons from PHAST - make sure you use the github version
	//
	chrom = maskedMaf.baseName.split('_')[0]
	"""
	${params.phast_path}./phastCons ${maskedMaf} ${autoNeutralModel} --rho ${params.rho} --target-coverage ${params.target_coverage} --expected-length ${params.expected_length} --no-post-probs --score --viterbi ${chrom}_phastcons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

process callNonAutosomalConservedElements {
	tag "Identifying conserved loci on non-autosomal ${fname}."

	publishDir params.outdir, mode: "copy", overwrite: true

	errorStrategy 'retry'
	maxRetries 3

	input:
	file(maskedMaf)

	output:
	file("*_unfilt.bed") 

	script:
	//
	// phastCons from PHAST - make sure you use the github version :/
	//
	chrom = maskedMaf.baseName.split('_')[0]
	"""
	${params.phast_path}./phastCons ${maskedMaf} ${baseDir}/${params.nonauto_neutral_model}${chrom}.mod --rho ${params.rho} --target-coverage ${params.target_coverage} --expected-length ${params.expected_length} --no-post-probs --score --viterbi ${chrom}_phastcons_unfilt.bed --seqname ${chrom} -i MAF
	"""
}

/*
* filter phastCons by log-odds score
*/

process filterByScore {
	tag "Filtering phastcons from ${fname} by log odds score"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastCons)

	output:
	file("*.bed")

	script:
	chrom = phastCons.baseName.split('_')[0]
	"""
	python ${baseDir}/bin/score_filter_phastcons.py ${phastCons} ${chrom}_phastcons_score_filt.bed ${params.min_decile}
	"""
}

/*
* filter for phastCons that are in syntenic
* regions based on with user-identified inputs
*/

process syntenyFilterPhastCons {
	tag "Synteny filtering phastCons from ${chrom}"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastCons)
	file(syntenyFilterFile)

	output:
	file("*phastcons_synteny_filt.bed")

	script:
	chrom  = phastCons.baseName.split('_')[0]
	"""
	bedtools intersect -a ${phastCons} -b ${syntenyFilterFile}  | bedtools sort -i - | bedtools merge -d 0 -i - > ${chrom}_phastcons_synteny_filt.bed
	"""
}

/*
* filter phastCons for elements that could confound the acceleration analysis
* e.g. duplicated regions, pseudogenes, etc.
*/

process phastConsDupFilter {
	tag "Filtering phastCons from ${fname} for repeats, duplicated regions, pseudogenes and selfChain"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(sizeFiltPhastCons)
	file(arFilterFile)

	output:
	file("*.bed")

	script:
	chrom = sizeFiltPhastCons.baseName.split('_')[0]
	"""
	bedtools subtract -a ${sizeFiltPhastCons} -b ${arFilterFile} | bedtools sort -i - | bedtools merge -d 0 -i - > ${chrom}_phastcons_dup_filt.bed
	"""
}


/*
* filter phastCons elements for minimum size of element
*/

process phastconsSizeFilter {
	tag "Filtering phastCons from ${fname} for size > 50bp"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastCons) 

	output:
	file("*.bed")

	script:
	chrom = file(phastCons).baseName.split('_')[0]
	"""
	awk '(\$3-\$2) >= 50' ${phastCons} > ${chrom}_phastcons_size_filt.bed
	"""
}


//* add ID to phastCons, necessary for extracting MAFs


process idToPhastcons {
	tag "Adding ID to phastCons from ${fname}"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(phastCons)

	output:
	file("*.bed") 

	script:
	chrom = phastCons.baseName.split('_')[0]
	"""
	awk -F'\t' -v OFS='\t' '{ \$(NF+1)=\$1"_"NR} 1' ${phastCons} > ${chrom}_phastCons_id.bed
	"""

}

/*
* run phyloP to identify accelerated elements, autosomes
*/

process accRegionsAutosomal {
	tag "Identifying accelerated regions on ${window}"

	errorStrategy 'retry'
	maxRetries 3

	input:
	tuple val(window), file(phastcons), file(species_maf) 

	output:
	stdout 

	script:
	//
	// phyloP from PHAST
	//
	chrom = phastcons.baseName.split('_')[0]
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -d ${params.random_seed} -g ${autoNeutralModel} ${species_maf}
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
	tuple val(window), file(phastcons), file(species_maf)

	output:
	stdout 

	script:
	//
	// phyloP from PHAST
	//
	chrom = phastcons.baseName.split('_')[0]
	"""
	${params.phast_path}./phyloP --features ${phastcons} --msa-format MAF --method LRT --mode ACC --subtree ${params.species_of_interest} -d ${params.random_seed} -g ${baseDir}/${params.nonauto_neutral_model}${chrom}.mod ${species_maf}
	"""

}

/*
* multiple test correction with Benjamini-Hochberg methodology
*/


process multipleTestCorrectionFiltering {
	tag "Multiple test correction and filtering"

	publishDir params.outdir, mode: "copy", overwrite: true

	input:
	file(scored_phastcons) 

	output:
	file("final_ARs_${params.random_seed}.bed")
	file("scored_phastCons_${params.random_seed}.txt")

	script:
	"""
	python ${baseDir}/bin/multiple_test_correction_filtering.py ${scored_phastcons} ${params.max_p} final_ARs_${params.random_seed}.bed scored_phastCons_${params.random_seed}.txt
	"""

}

speciesFile = file(params.species)
chromBedPath = file(params.chrom_bed_path)
autoNeutralModel = file(params.auto_neutral_model)
arFilterFile = file(params.ar_filters_path)

// combine synteny files into one file for easier filtering
// Channel.fromPath( "${params.synteny_filter_path}*.bed", checkIfExists: true ).collect().flatten().collectFile(name: 'synteny_filters.bed', storeDir: params.outdir, sort: false)
syntenyFilterFile = file(params.synteny_filter_file)

workflow {

	// set up channel per chromosome
	Channel.fromPath( "${params.maf_path}chr*.maf.gz" ).set{ initMaf }

	// filter MAF to only include species of interest
	extractSpecies(speciesFile, initMaf).set{ speciesMaf }

	// mask species of interest 
	maskSpeciesOfInterest(speciesMaf, chromBedPath).set{ maskedMaf }

	// run phastCons to get conserved elements (not including masked sequence)
	// this is done separately for autosomal vs. nonautosomal because of the 
	// different neutral models
	callAutosomalConservedElements(maskedMaf.filter( ~/.*chr\d.*/ ), autoNeutralModel).set{ autoPhastConsElements }
	callNonAutosomalConservedElements(maskedMaf.filter( ~/.*chr[X Y M].*/ )).set{ nonAutoPhastConsElements }

	// combine autosomal and non-autosomal phastCons elements
	phastConsElements = autoPhastConsElements.concat(nonAutoPhastConsElements)

	// filter phastCons elements by log-odds score to get the most conserved elements
	filterByScore(phastConsElements).set{ scoreFiltPhastCons } 

	// filter phastCons elements to keep only those that are syntenic
	// in selected species, traditionally rhesus, dog, mouse relative to human
	// based on concatenated, sorted BED version of net files downloaded from UCSC Table Browser
	syntenyFilterPhastCons(scoreFiltPhastCons, syntenyFilterFile).set{ syntFiltPhastCons }

	// filter phastCons based on ENCODE blocklist
	phastConsDupFilter(syntFiltPhastCons, arFilterFile).set{ phastConsDupFiltered }

	// filter by minimum element size
	phastconsSizeFilter(phastConsDupFiltered).set{ phastConsSizeFilt }

	// add ID (needed to extract MAFs for phyloP)
	idToPhastcons(phastConsSizeFilt).set{ phastConsFiltered }

	// run phyloP to identify accelerated elements, autosomes
	// need to coordinate phastcons with relevant MAF based on chromosome
	phastConsSplitAuto = phastConsFiltered.flatten().filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.baseName.split(/_phastCons/)[0], file)
	}

	phastConsSplitNonAuto = phastConsFiltered.flatten().filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.baseName.split(/_phastCons/)[0], file)
	}

	speciesMafChromAuto = speciesMaf.filter( ~/.*chr\d.*/ ).map {
	file -> tuple(file.baseName.split(/.maf/)[0], file)
	}

	speciesMafChromNonAuto = speciesMaf.filter( ~/.*chr[X Y M].*/ ).map {
	file -> tuple(file.baseName.split(/.maf/)[0], file)
	}

	phastConsSplitAutoMaf = phastConsSplitAuto.combine(speciesMafChromAuto, by: 0)

	phastConsSplitNonAutoMaf = phastConsSplitNonAuto.combine(speciesMafChromNonAuto, by: 0)

	// phyloP for acceleration, then multiple test correction of p-values 
	accRegionsAutosomal(phastConsSplitAutoMaf).set{ autosomalARs }
	accRegionsNonAutosomal(phastConsSplitNonAutoMaf).set{ nonAutosomalARs }
	acceleratedRegions = autosomalARs.merge(nonAutosomalARs).collect().flatten()//.collectFile(name: 'phyloP_results.txt')
	acceleratedRegions.collectFile(name: "phyloP_results.txt", storeDir: params.outdir, sort:false).set{ acceleratedRegionsFile }

	multipleTestCorrectionFiltering(acceleratedRegionsFile)

}



