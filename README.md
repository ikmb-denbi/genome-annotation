![](images/ikmb_bfx_logo.png) ![](images/deNBI_Logo_rgb.jpg)
# Genome Annotation 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](http://singularity.lbl.gov)

This pipeline can be used to annotate a genome *de novo*. 

### Pipeline main steps  

1. Hints file is generated from all available evidences (proteins/EST and/or RNA-seq reads). 

2. Gene models are predicted using Augustus with the hints file as extrinsic evidence. 

3. Gene models are functionally annotated using Blast and InterProScan. 

The pipeline runs these three steps by default, but each step can also be skipped by command line parameters. The minimum requirement is at least one type of evidence. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. Most of the required programs are installed via conda during the execution, making installation trivial and results highly reproducible.

### Documentation 

Documentation about the pipeline can be found in the `docs/` directory or under the links below:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits 

This pipeline was written by Dr. Montserrat Torres ([m.torres](https://git.ikmb.uni-kiel.de/m.torres)) at [IKMB](http://www.ikmb.uni-kiel.de).
