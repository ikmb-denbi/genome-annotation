# Genome Annotation - Installation

## 1. Install Miniconda2

In the IKMB clusters:  

`module load miniconda2` 

Otherwise, install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/)

## 2. Create and load the Hints environment 

The file `environment.yml` in the main directory contains instructions to install all the necessary programs to run the pipeline in a conda environment called "Hints". 

To create this environments type the following command: 

`conda env create -f environment.yml` 

To activate the environment type: 

`source activate Hints` 

All the necessary programs (except GenomeThreader, see below) are now installed so you can run the NF-hints pipeline. 

## 3. Install GenomeThreader 

If you want to use GenomeThreader (in addition to Exonerate) to map the proteins file to your genome and create additional hints, you have to install it yourself. 

[GenomeThreader download](http://genomethreader.org/download.html) 

Make sure the program is in your path when you run the pipeline. 

By default, the pipeline will try to run GenomeThreader. If it is not installed you will get an error message. To run the pipeline without GenomeThreader add the argument: 

`--gth false` 


