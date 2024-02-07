# AcceleratedRegionsNF
Pipeline to identify accelerated regions given a multiple alignment of genomes.

Updated to Nextflow DSL 2 in Feb. 2024. Please report any bugs you find by opening an issue. Thanks!

Please refer to [Pollard et al. 2006](https://www.nature.com/articles/nature05113), [Lindblad-Toh et al. 2011](https://www.nature.com/articles/nature10530) and [Hubisz & Pollard 2014](https://www.sciencedirect.com/science/article/pii/S0959437X14000781) for more information on accelerated regions and methods behind this pipeline. For more information on Nextflow please [refer to its documentation](https://www.nextflow.io/docs/latest/index.html).

**How to call accelerated regions with this pipeline:**

1.) clone the repo

2.) install the required packages (see [wiki](https://github.com/keoughkath/AcceleratedRegionsNF/wiki/Required-packages))

3.) make a copy of and adapt the project file (example: `ARs_config_20way_example_human.yml`) to fit your project goals
* <ins>species</ins>: a newline-delimited list of species in a text file for inclusion in the analysis, this can be a subset of or all the species in your multiple alignment
* <ins>maf_path</ins>: the path to your multiple alignment files which must be in [MAF format](https://genome.ucsc.edu/FAQ/FAQformat.html#format5)
* <ins>phast_path</ins>: path to your installation of PHAST, this is specified because we recommend using [the Github version of PHAST](https://github.com/CshlSiepelLab/phast) and this way you don't have to alter your $PATH
* <ins>outdir</ins>: where you want the outputted files to be saved
* <ins>init_tree</ins>: a rooted, bifurcating, Newick-formatted tree
* <ins>species_of_interest</ins>: the species in which you want to identify lineage-specific accelerated elements
* <ins>chrom_bed_path</ins>: BED files specifying the size of each chromosome analyzed in the reference frame of the MAF reference sequence, these should be named as e.g. "chrX.bed"
* <ins>synteny_filter_file</ins>: BED-formatted (not gzipped) file syntenic regions that will be intersected with the phastCons
* <ins>ar_filters_path</ins>: path to a gzipped BED file containing anything you want to exclude from the accelerated region analysis, e.g. repetitive elements, blocklists, etc.
* <ins>max_p</ins>: Maximum Benjamini-Hochberg-adjust p-value you would like to consider for accelerated regions (the scored phastCons elements are also returned so you can change your mind about this value posthoc)
* <ins>random_seed</ins>: random seed used for the phyloP step for reproducibility of results
* <ins>target_coverage</ins>: phastCons parameter (aka gamma) that indicates the percent of the query genome expected to be conserved
* <ins>expected_length</ins>: phastCons parameter (aka omega) that indicates the expected length of conserved elements
* <ins>rho</ins>: phastCons scaling parameter indicating how to scale the neutral tree to obtain the conserved tree (see [Siepel et al. 2005](https://pubmed.ncbi.nlm.nih.gov/16024819/) for more details)
* <ins>auto_neutral_model</ins>: path to the neutral model for autosomes (not ending in "X", "Y", or "M", so this may need to be adjusted for some species) in ".mod" format
* <ins>nonauto_neutral_model</ins>: path to a directory containing the neutral models for non-autosomes (ending in "X", "Y", or "M", so this may need to be adjusted for some species) in ".mod" format, these should be named e.g. "chrX.mod"
* <ins>min_decile</ins>: length-normalized log odds score threshold for phastCons elements, a higher decile indicates higher conservation (in general) see [Siepel et al. 2005](https://pubmed.ncbi.nlm.nih.gov/16024819/) and [phastCons documentation](http://compgen.cshl.edu/phast/help-pages/phastCons.txt) for more details

4.) adjust the nextflow config file to match your operating environment (we have provided a sample config for an SGE system)

5.) run the pipeline (sample command: `nextflow run call_ARs.nf -w "hars_workdir/" -profile local -params-file ARs_config_20way_example_human.yml`)

The example files indicated are for the UCSC 20-way primate alignment, which runs on a 2020 MacBook Pro in a few minutes, and is good for testing.

**Outputs**

If the pipeline successfully runs, it will output two files, final_ARs_(random_seed).bed and scored_phastCons_(random_seed).txt, which are a BED file of the accelerated regions and a tab-separated file of the filtered, acceleration-scored phastCons which can be useful in statistical analyses as null models or if you want to reselect HARs at a different FDR. A bunch of files are also outputted by the intermediate steps, which you can turn off if desired. 
