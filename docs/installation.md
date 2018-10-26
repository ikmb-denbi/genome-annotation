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

To run the pipeline you first have to check out the code to a location in your system (i.e. $HOME/git/). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
## Installing all other software 

This pipeline uses a lot of different bioinformatics software (you can see a full list at the end of this document). How you proceed to install these programs will depend on what system/cluster you are using: 

### 1. Working in the IKMB RZ cluster 

In the **IKMB RZ cluster**, all these programs are available as modules and will be loaded automatically as the pipeline runs. You don't need to do anything else, just make sure that you run the pipeline using the parameter `-profile standard`. This is anyway the default profile, so you don't even need to especify it, only don't use any other. 

### 2. Not in the IKMB? Use Bioconda 

Most of the required programs are available as [bioconda packages](https://bioconda.github.io/recipes.html) for easy installation. All you need to do is install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/) 

In this case, when you run the pipeline add `-profile conda` to your command. 

Nextflow will install the environments with the necessary packages as it runs. 

#### Missing programs: 

Some of the required programs are not available as conda packages yet, so you will have to install them yourself. Some parts of the pipeline will run successfully anyway, but you need to turn off the ones that don't. 

1. **GenomeThreader:** it is used to map the protein evidences to the genome and create hints. It is not necessary to run the pipeline and, by default, protein evidence will be mapped using [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate). However, my experience shows that the two programs complement each other: some models that were missed by Exonerate are found by GenomeThreader, and viceversa. You can download GenomeThreader from [here](http://genomethreader.org/download.html) and follow the installation instructions. If you want to run the pipeline without this program you must use `--gth false` on your command call.

2. **Annie:** it is used to transfer the functional annotations to your gff3 file. You can download it from [here](http://genomeannotation.github.io/annie/) and extract it in your system. Make sure the executable file `annie.py` is in your path. If you don't want to perform functional annotation use `--funAnnot false`. 

3. **InterProScan:** it is also necessary during the functional annotation step. It it is used to scan the predicted protein sequences for known protein signatures (functional domains, GO terms, etc) searching in different public databases. You can download it from [here](https://www.ebi.ac.uk/interpro/download.html) and extract it in your system. Make sure the executable `interproscan.sh` is in your path. If you want to run the pipeline without this program use `--funAnnot false`. 

4. **Bioruby:** to run the functional annotation step you also need to have the bioruby library installed. You can download and install ruby if you haven't done so yet from [here](https://www.ruby-lang.org/en/). Then install bioruby using the RubyGems tool: 

`gem install bio` 


*coming soon:* [Singularity]() image with these programs installed. 

#### Configuration for other clusters

### 3. Install all programs yourself 

Here is a list of all the programs necessary to run the complete genome-annotation pipeline (--prots proteins.fa --ESTs ESTs.fa --reads '*_R{1,2}.fastq' --gth true --RM true --trinity true --augustus true --funAnnot true). The pipeline has been tested successfully  with the versions that are here described. 

1. Blast+ v2.2.30 

2. Exonerate v2.2.0 

3. Bioperl 

4. GenomeThreader v1.7.0 

5. Genometools v1.5.6 

6. RepeatMasker v.4.0.6 

7. Trimgalore v0.4.4 

8. Hisat2 v2.1.0 

9. Samtools v1.5 

10. Augustus v3.2.1 

11. Trinity v2.5.1 

12. Bowtie v2.2.3 

13. Ruby v2.2.2 + bioruby library 

14. Annie  

15. Interproscan 











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
