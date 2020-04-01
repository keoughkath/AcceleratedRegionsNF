# AcceleratedRegionsNF
Pipeline to identify accelerated regions given a multiple alignment of genomes.

This pipeline is new and likely has many bugs. Please report any you find by opening an issue. Thanks!

Please refer to [Pollard et al. 2006](https://www.nature.com/articles/nature05113), [Lindblad-Toh et al. 2011](https://www.nature.com/articles/nature10530) and [Hubisz & Pollard 2014](https://www.sciencedirect.com/science/article/pii/S0959437X14000781) for more information on accelerated regions and methods behind this pipeline. For more information on Nextflow please [refer to its documentation](https://www.nextflow.io/docs/latest/index.html).

How to call accelerated regions with this pipeline:

1.) clone the repo

2.) install the Conda environment contained in the repo (ARs_conda.yml)

3.) make a copy of and adapt the project.yml file to fit your project goals

4.) adjust the nextflow config file to match your operating environment (optionally add an additional config file to match a parallel system such as SGE)

5.) run the pipeline (sample command: `nextflow run call_ARs.nf -w "output_dir" -profile local -params-file project.yml`)
