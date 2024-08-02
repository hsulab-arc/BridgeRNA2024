# BridgeRNA2024

This code repository accompanies the paper ["Bridge RNAs direct programmable recombination of target and donor DNA"](https://www.nature.com/articles/s41586-024-07552-4) by Durrant & Perry et al, 2024.

It includes snakemake pipelines and code for various key analyses described in the paper.

# Contents
There are 5 separate pipelines in this code repository. Pipelines are meant to be run through the `brna2024` command.
Each pipeline includes a `*.Snakefile` and `*.config.yml` file. There is also a `*.py` file that contains the
python functions for the pipeline. The `*.py` file is imported by the `*.Snakefile` and is used to implement the
snakemake rules. There are also accompanying `*.R` files that contain R scripts that implement some of the rules.

### `rnaseq`
  * This pipeline was used to process the small RNA-seq data from WT donor plasmids of the IS621 and related orthologs 
  and aligns them to the donor plasmids to detect the expression of the bridge RNAs.
  * Related files
    * `rnaseq.py`
    * `snakemake/rnaseq.Snakefile`
    * `snakemake/rnaseq.config.yml`
### `structure`
  * This pipeline was used to detect the secondary structure of bridge RNAs in the flanks of the bridge 
  recombinase CDS by the alignment of many homologous sequences.
  * Related files
    * `structure.py`
    * `snakemake/structure.Snakefile`
    * `snakemake/structure.config.yml`
### `targetscreen`
  * This pipeline was used to process the amplicon sequencing data from the target screen experiments.
  * Related files
    * `targetscreen.py`
    * `snakemake/targetscreen.Snakefile`
    * `snakemake/targetscreen.config.yml`
### `donorscreen`
  * This pipeline was used to process the amplicon sequencing data from the donor screen experiments.
  * Related files
    * `donorscreen.py`
    * `snakemake/donorscreen.Snakefile`
    * `snakemake/donorscreen.config.yml`
### `nanoporepipeline`
  * This pipeline was used to process the nanopore sequencing data from the genomic insertion experiments.
  * Related files
    * `nanoporepipeline.py`
    * `snakemake/nanoporepipeline.Snakefile`
    * `snakemake/nanoporepipeline.config.yml`

# Note
This repository is for reference only. It is not guaranteed to work on your system, and will likely require 
modifications to run properly. The Dockerfile can be useful for installing some of the dependencies. Once working,
it is intended to be used as follows:
```
brna2024 <pipeline name> --threads <number of threads> <workdir>
```
or if you want to run the pipeline through docker:
```
source setup.sh;
brna2024_docker_run <pipeline name> <number of threads> <workdir>
```

# Citation
Please cite [our paper](https://www.nature.com/articles/s41586-024-07552-4) if you use any aspect of this 
code repository in your work.
