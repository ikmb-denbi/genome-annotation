![](images/deNBI_logo.jpg) ![](images/ikmb_bfx_logo.png) 

# ESGA - Genome Annotation 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](http://singularity.lbl.gov)

This pipeline can be used to annotate a genome *de novo*. 

### Pipeline main steps  

1. Hints file is generated from all available evidences (proteins, EST and/or RNA-seq reads). 

2. Gene models are predicted using Augustus with the hints file as extrinsic evidence. 

The pipeline runs these three steps by default, but each step can also be skipped by command line parameters. The minimum requirements are a genome file and at least one type of evidence. 

### Test data

A simple test data set can be downloaded [here](https://drive.google.com/open?id=1VFqLnRJiuj5Vhj2KCOdY58jwxZKkkMVU) 
 
### Documentation 

Documentation about the pipeline can be found in the `docs/` directory or under the links below:

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)
6. [What's next](docs/whatsnext.md)

### Pipeline Scheme 

![](images/genome-annotation_dag_mod.svg) 

### Credits 

This pipeline was written by Dr. Montserrat Torres ([MontseTor](https://github.com/MontseTor)) and Dr. Marc Höppner ([marchoeppner](https://github.com/marchoeppner)) at [IKMB](http://www.ikmb.uni-kiel.de).
