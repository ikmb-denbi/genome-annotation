!!! THIS PROJECT HAS BEEN RETIRED. PLEASE WORK WITH ITS DROP-IN REPLACEMENT: https://github.com/ikmb/esga

![](images/deNBI_logo.jpg) ![](images/ikmb_bfx_logo.png) 

# ESGA - Genome Annotation 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](http://singularity.lbl.gov)

This pipeline can be used to generate "hints" from aligned sequence evidence to annotate a genome *de novo*. 

### Pipeline main steps  

1. Hints file is generated from all available evidences (proteins, EST and/or RNA-seq reads). 

2. Gene models are predicted using Augustus with the hints file as extrinsic evidence (optional). 

The minimum requirements are a genome file and at least one type of evidence. 

### Test data

A simple test data set can be downloaded [here](https://drive.google.com/open?id=1VFqLnRJiuj5Vhj2KCOdY58jwxZKkkMVU) 
 
### Documentation 

Documentation about the pipeline can be found in the `docs/` directory or under the links below:

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Recommendations](docs/recommendations.md)
3. [Installation and configuration](docs/installation.md)
4. [Running the pipeline](docs/usage.md)
5. [Output](docs/output.md)
6. [Troubleshooting](docs/troubleshooting.md)
7. [What's next](docs/whatsnext.md)

### Pipeline Scheme 

![](images/Pipeline_dag.svg) 

### Credits 

This pipeline was written by Dr. Montserrat Torres ([MontseTor](https://github.com/MontseTor)) and Dr. Marc Höppner ([marchoeppner](https://github.com/marchoeppner)) at [IKMB](http://www.ikmb.uni-kiel.de).
The authors gratefully acknowledge inspiration, fruitful discussions and a few useful code snippets from the [nf-core](https://www.nf-co.re).
