# Installation and configuration 

## Installing Nextflow 

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
Therefore the first thing to do is install Nextflow. 

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

You can check out the code to a location on your system (i.e. $HOME/git/). This is recommended as you will have to make some minor changes so the pieline knows how to run on your system (see below). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
Launching the pipeline then works as follows:

```bash
nextflow run $HOME/git/genome-annotation/main.nf <run options>
```

## Creating a config file for your cluster

This pipeline uses at minimum two configuration files. The file `conf/base.config` contains information about the resource requirements of the individual stages of the pipeline. Normally, you do not have to touch this.

In addition, you will need a config file that specifies which resource manager your cluster uses. An example for a Slurm cluster which uses Singularity (see below) is included as [conf/slurm_ikmba.config](../conf/slurm_ikmba.config). 
Detailed instructions about resource managers and available options can be found [here](https://www.nextflow.io/docs/latest/executor.html).

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

This pipeline uses a lot of different bioinformatics software - you can find a full list with versions in the included file [environment.yml](../environment.yml).
You won't have to install of of these tools, but can instead use one of the two options below:

### A. (Bio-) Conda

Most of the required programs are available as [bioconda packages](https://bioconda.github.io/recipes.html) for easy installation. All you need to do is install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/) 

In this case, when you run the pipeline add `-profile conda` (if in the IKMB) or `-profile conda_custom` (see below) to your command. 

However, please note that you will not be able to run the GenomeThreader part of the pipeline unless you make sure to provide this package some other way. 

Nextflow will install the environments with the necessary packages during pipeline start-up. 

### B. Singularity

The preferred way of provisioning the software is through [Singularity](https://github.com/sylabs/singularity). If Singularity is not available on your cluster, please ask your admins to install it. 

To enable use of singularity, simply add the following to your custom config file (see below):

```bash
singularity {
	enabled = true
}
```

Depending on your cluster and configuration of singularity, you may also have to provide some additional run options. A typical example would be that your data is stored on a network-mounted drive, which is not automatically detected by singularity. In this case, you can do:

```bash
singularity {
	enabled=true
	runOptions="-B /path/to/network/drive"
}
```
