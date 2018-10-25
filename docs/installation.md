# Genome Annotation - Installation and configuration

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


# Genome Annotation - Local Configuration

If running the pipeline in a local environment, we highly recommend using either Docker or Singularity.

## Docker
Docker is a great way to run NF-hints, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required. The NF-hints profile comes with a configuration profile for docker, making it very easy to use. This also comes with the required presets to use the AWS iGenomes resource, meaning that if using common reference genomes you just specify the reference ID and it will be autaomtically downloaded from AWS S3.

First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run NF-hints -profile docker --reads '<path to your reads>'
```

Nextflow will recognise `NF-hints` and download the pipeline from GitHub. The `-profile docker` configuration lists the [NF-hints](https://hub.docker.com/r/NF-hints/) image that we have created and is hosted at dockerhub, and this is downloaded.

For more information about how to work with reference genomes, see [`docs/configuration/reference_genomes.md`](docs/configuration/reference_genomes.md).

### Pipeline versions
The public docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r v1.3`. This uses pipeline code and docker image from this tagged version.


## Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity 'docker://NF-hints'`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name NF-hints.img docker://NF-hints
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/NF-hints -with-singularity /path/to/NF-hints.img
```
