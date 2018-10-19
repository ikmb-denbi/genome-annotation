# Genome Annotation 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](http://singularity.lbl.gov)

This pipeline can be used to annotate a genome de novo. 

### Introduction 

First, a Hints file is generated from all available evidences (proteins/EST and/or RNA-seq reads). 

Second, Augustus is used to predict gene models in the genome using the Hints file as extrinsic evidence. 

Finally, functional annotation can also be run. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

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

This pipeline was written by MTorres ([m.torres](https://git.ikmb.uni-kiel.de/m.torres)) at [IKMB](http://www.ikmb.uni-kiel.de).
