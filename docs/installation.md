# Installation and configuration 

## Compute infrastructure

This pipeline is designed to run on a distributed compute system, such as a traditional HPC cluster system. 
We have tested the pipeline on two Slurm clusters, with node configurations of 16cores/128GB Ram and 20cores/256GB Ram, respectively. 

While smaller nodes will probably work, it may require some tweaking on your end. Most importantly, if you plan on using the transcriptome 
assembly branch of the pipeline, available memory may become limiting (however, 128GB Ram should be fine for typical datasets; 256GB are perhaps 
necessary if you plan on using a larger sample size). 

## Installing Nextflow 

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across 
distributed compute systems in a very portable manner. Therefore the first thing to do is to install Nextflow. 

Nextflow runs on most systems (Linux, Mac OSX etc). It can be installed by running the following commands:

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

**You need NextFlow version >= 19.03.0 to run this pipeline.** 

## Obtaining the code 

You can check out the code to a location on your system (i.e. $HOME/git/). This is recommended as you will have to make some minor changes 
so the pieline knows how to run on your system (see below). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
Launching the pipeline then works as follows:

```bash
nextflow run $HOME/git/genome-annotation/main.nf <run options>
```

## Creating a config file for your cluster

This pipeline uses at minimum two configuration files. The file [conf/base.config](../conf/base.config) contains information about the resource requirements 
of the individual stages of the pipeline. Normally, you do not have to touch this.

In addition, you will need a config file that specifies which resource manager your cluster uses. An example for a Slurm cluster which uses 
Singularity (see below) is included as [conf/slurm_ikmba.config](../conf/slurm_ikmba.config). Detailed instructions about resource managers and 
available options can be found [here](https://www.nextflow.io/docs/latest/executor.html).

You can now create a new execution profile by editing the file [nextflow.config](../nextflow.config):

```bash
profiles {

	...
	your_cluster {
		includeConfig 'conf/base.config'
		includeConfig 'conf/your_cluster.config'
	}
}
```

## Installing all other software 

This pipeline uses a lot of different bioinformatics software - you can find a full list with versions in the included 
file [environment.yml](../environment.yml). You won't have to install any of these tools, but can instead use one of the two options below:

### A. (Bio-) Conda

All of the required programs are available as [bioconda packages](https://bioconda.github.io/recipes.html) for easy installation. 
All you need to do is install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/) 

Your custom config file should then contain the following line:

```bash
process {
	conda = "${baseDir}/environment.yml"
}
```

Nextflow will install the environment with the necessary packages during pipeline start-up. However, please be advised that this takes a little while and 
needs to be done every time you run a new project. Using Singularity is therefore highly recommended (see below). Failing that, you can also install the conda environment into some other location (default is usually $HOME/.conda) and just make sure to activate it any time you want to start the pipeline. Your config file would then container neither a `conda` nor a `singularity` statement. 

### B. Singularity

The preferred way of provisioning the software is through [Singularity](https://github.com/sylabs/singularity). If Singularity is not available on 
your cluster, please ask your admins to install it. 

To enable use of singularity, simply add the following to your custom config file (see below):

```bash
singularity {
	enabled = true
}
```

Depending on your cluster and configuration of singularity, you may also have to provide some additional run options. 
A typical example would be that your data is stored on a network-mounted drive, which is not automatically detected by singularity. In this case, you can do:

```bash
singularity {
	enabled=true
	runOptions="-B /path/to/network/drive"
}
```
