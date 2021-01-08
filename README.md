# AcceleratedRegionsNF
Pipeline to identify accelerated regions given a multiple alignment of genomes.

This pipeline is new and likely has many bugs. Please report any you find by opening an issue. Thanks!

Please refer to [Pollard et al. 2006](https://www.nature.com/articles/nature05113), [Lindblad-Toh et al. 2011](https://www.nature.com/articles/nature10530) and [Hubisz & Pollard 2014](https://www.sciencedirect.com/science/article/pii/S0959437X14000781) for more information on accelerated regions and methods behind this pipeline. For more information on Nextflow please [refer to its documentation](https://www.nextflow.io/docs/latest/index.html).

Changes as of 8 January 2021:
* the neutral model is now a parameter you add, so you can use a neutral model defined by multiple methods, e.g. 4D sites or ancestral repeats
* this pipeline will now work on MAFs chunked e.g. into 10Mb blocks, as the Zoonomia files are formatted, hopefully saving on memory during compute

How to call accelerated regions with this pipeline:

1.) clone the repo

2.) install the Conda environment contained in the repo (ARs_conda.yml)

3.) make a copy of and adapt the project file (example: zoonomia_example_human.yml) to fit your project goals

4.) adjust the nextflow config file to match your operating environment (optionally add an additional config file to match a parallel system such as SGE)

5.) run the pipeline (sample command: `nextflow run call_ARs_zoo.nf -w "zoonomia/zoo_hars_full/" -profile sge -params-file human_all_species_zoo_full_run2.yml`)
