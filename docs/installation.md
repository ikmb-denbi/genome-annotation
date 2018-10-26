# Installation and configuration 

## Installing Nextflow 

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
Therefore the first thing to do is install Nextflow. 

If you are working in the **IKMB RZ cluster** you can simply load the following modules:

`module load IKMB Java/1.8.0 Nextflow/0.32.0` 

Otherwise, you can easily install it yourself. Nextflow runs on most systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

**You need NextFlow version >= 0.32 to run this pipeline.** 

## Cloning the genome-annotation repository 

To run the pipeline you first have to check out the code to a location where you have read and write permissions (i.e. $HOME/git/). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
## Installing all other software 

This pipeline uses a lot of different bioinformatics software (you can see a full list at the end of this document). How you proceed to install these programs will depend on what system/cluster you are using: 

### 1. Working in the IKMB RZ cluster 

In the **IKMB RZ cluster**, all these programs are available as modules and will be loaded automatically as the pipeline runs. You don't need to do anything else, just make sure that you run the pipeline using the parameter `-profile standard`. This is anyway the default profile, so you don't even need to especify it, only don't use any other. 

### 2. Not in the IKMB? Use Bioconda 

Most of the required programs are available as [bioconda packages](https://bioconda.github.io/recipes.html) for easy installation. 









In the IKMB clusters:  

`module load miniconda2` 

Otherwise, install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/)

## 2. Create and load the Hints environment 

The file `environment.yml` in the main directory contains instructions to install all the necessary programs to run the pipeline in a conda environment called “Hints”. 

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
nextflow run NF-hints -profile docker --reads ‘<path to your reads>‘
```

Nextflow will recognise `NF-hints` and download the pipeline from GitHub. The `-profile docker` configuration lists the [NF-hints](https://hub.docker.com/r/NF-hints/) image that we have created and is hosted at dockerhub, and this is downloaded.

For more information about how to work with reference genomes, see [`docs/configuration/reference_genomes.md`](docs/configuration/reference_genomes.md).

### Pipeline versions
The public docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r v1.3`. This uses pipeline code and docker image from this tagged version.


## Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity ‘docker://NF-hints’`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you’ll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name NF-hints.img docker://NF-hints
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/NF-hints -with-singularity /path/to/NF-hints.img
```

# genome-annotation: Configuration for other clusters

It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the pipeline knows how to work with your cluster.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config` (see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more).

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile - see [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

## Cluster Environment
By default, pipeline uses the `local` Nextflow executor - in other words, all jobs are run in the login session. If you’re using a simple server, this may be fine. If you’re using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```nextflow
process {
  executor = ‘YOUR_SYSTEM_TYPE’
}
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```nextflow
process {
  executor = ‘SLURM’
  clusterOptions = ‘-A myproject’
}
```


## Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system. If possible, we _highly_ recommend using either Docker or Singularity.

### Docker
Docker is a great way to run NF-hints, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required.

First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, simply run the analysis pipeline:
```bash
nextflow run NF-hints -profile docker --reads ‘<path to your reads>‘
```

Nextflow will recognise `NF-hints` and download the pipeline from GitHub. The `-profile docker` configuration lists the [NF-hints](https://hub.docker.com/r/NF-hints/) image that we have created and is hosted at dockerhub, and this is downloaded.

The public docker images are tagged with the same version numbers as the code, which you can use to ensure reproducibility. When running the pipeline, specify the pipeline version with `-r`, for example `-r v1.3`. This uses pipeline code and docker image from this tagged version.

To add docker support to your own config file (instead of using the `docker` profile, which runs locally), add the following:

```nextflow
docker {
  enabled = true
}
process {
  container = wf_container
}
```

The variable `wf_container` is defined dynamically and automatically specifies the image tag if Nextflow is running with `-r`.

A test suite for docker comes with the pipeline, and can be run by moving to the [`tests` directory](https://github.com/NF-hints/tree/master/tests) and running `./run_test.sh`. This will download a small yeast genome and some data, and attempt to run the pipeline through docker on that small dataset. This is automatically run using [Travis](https://travis-ci.org/NF-hints/) whenever changes are made to the pipeline.

### Singularity image
Many HPC environments are not able to run Docker due to security issues. [Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker. Even better, it can use create images directly from dockerhub.

To use the singularity image for a single run, use `-with-singularity ‘docker://NF-hints’`. This will download the docker container from dockerhub and create a singularity image for you dynamically.

To specify singularity usage in your pipeline config file, add the following:

```nextflow
singularity {
  enabled = true
}
process {
  container = “docker://$wf_container”
}
```

The variable `wf_container` is defined dynamically and automatically specifies the image tag if Nextflow is running with `-r`.

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you’ll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name NF-hints.img docker://NF-hints
```

Then transfer this file and run the pipeline with this path:

```bash
nextflow run /path/to/NF-hints -with-singularity /path/to/NF-hints.img
```


### Manual Installation
As a last resort, you may need to install the required software manually. We recommend using [Bioconda](https://bioconda.github.io/) to do this. The following instructions are an example only and will not be updated with the pipeline.

#### 1) Install miniconda in your home directory
``` bash
cd
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2) Add the bioconda conda channel (and others)
```bash
conda config --add channels anaconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels salilab
```

#### 3) Create a conda environment, with all necessary packages:
```bash
conda create --name NF-hints_py2.7 python=2.7
source activate NF-hints_py2.7
conda install --yes \
    fastqc \
    multiqc
```
_(Feel free to adjust versions as required.)_

##### 4) Usage
Once created, the conda environment can be activated before running the pipeline and deactivated afterwards:

```bash
source activate NF-hints_py2.7
# run pipeline
source deactivate
```

## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS=’-Xms1g -Xmx4g’
```



## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run NF-hints --reads ‘*_R{1,2}.fastq.gz’ -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you’re running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull NF-hints
```

### Reproducibility
It’s a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you’ll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [NF-hints releases page](https://github.com/NF-hints/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you’ll know what you used when you look back in the future.
