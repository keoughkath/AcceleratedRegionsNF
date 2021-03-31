# AcceleratedRegionsNF
Pipeline to identify accelerated regions given a multiple alignment of genomes.

This pipeline is new and likely has many bugs. Please report any you find by opening an issue. Thanks!

Please refer to [Pollard et al. 2006](https://www.nature.com/articles/nature05113), [Lindblad-Toh et al. 2011](https://www.nature.com/articles/nature10530) and [Hubisz & Pollard 2014](https://www.sciencedirect.com/science/article/pii/S0959437X14000781) for more information on accelerated regions and methods behind this pipeline. For more information on Nextflow please [refer to its documentation](https://www.nextflow.io/docs/latest/index.html).

**How to call accelerated regions with this pipeline:**

1.) clone the repo

2.) install the Conda environment contained in the repo (ARs_conda.yml)

3.) make a copy of and adapt the project file (example: zoonomia_example_human.yml) to fit your project goals
* species: a newline-delimited list of species in a text file for inclusion in the analysis, this can be a subset of or all the species in your multiple alignment
* maf_path: the path to your multiple alignment files which must be in [MAF format](https://genome.ucsc.edu/FAQ/FAQformat.html#format5)
* phast_path: path to your installation of PHAST, this is specified because we recommend using the Github version of PHAST and this way you don't have to alter your $PATH
* outdir: where you want the outputted files to be saved
* init_tree: a rooted, bifurcating tree
* species_of_interest: the species in which you want to identify lineage-specific accelerated elements
* chrom_bed_path: BED files specifying the size of each chromosome analyzed in the reference frame of the MAF reference sequence, these should be named as e.g. "chrX.bed"
* synteny_filter_path: path to the directory containing BED-formatted (not gzipped) syntenic regions that will be intersected with the phastCons (1 or more files should be in this directory)
* ar_filters_path: path to a gzipped BED file containing anything you want to exclude from the accelerated region analysis, e.g. repetitive elements, blocklists, etc.
* max_p: Maximum Benjamini-Hochberg-adjust p-value you would like to consider for accelerated regions (the scored phastCons elements are also returned so you can change your mind about this value posthoc)
* random_seed: random seed used for the phyloP step for reproducibility of results
* target_coverage: phastCons parameter (aka gamma) that indicates the percent of the query genome expected to be conserved
* expected_length: phastCons parameter (aka omega) that indicates the expected length of conserved elements
* rho: phastCons scaling parameter indicating how to scale the neutral tree to obtain the conserved tree (see [Siepel et al. 2005](https://pubmed.ncbi.nlm.nih.gov/16024819/) for more details)
* auto_neutral_model: path to the neutral model for autosomes (not ending in "X", "Y", or "M", so this may need to be adjusted for some species) in ".mod" format
* nonauto_neutral_model: path to the neutral model for non-autosomes (ending in "X", "Y", or "M", so this may need to be adjusted for some species) in ".mod" format
* min_decile: length-normalized log odds score threshold for phastCons elements, a higher decile indicates higher conservation (in general) see [Siepel et al. 2005](https://pubmed.ncbi.nlm.nih.gov/16024819/) and [phastCons documentation](http://compgen.cshl.edu/phast/help-pages/phastCons.txt) for more details

4.) adjust the nextflow config file to match your operating environment (we have provided a sample config for an SGE system)
* change the path to "process.conda" so it matches the location of your installation of the `ARs_conda_env`
* adjust the paths for `-e` and `-o`

5.) run the pipeline (sample command: `nextflow run call_ARs_zoo.nf -w "hars_workdir/" -profile sge -params-file zoonomia_example_human.yml`)
