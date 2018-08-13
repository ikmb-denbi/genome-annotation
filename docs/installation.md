# NF-hints Installation

To start using the NF-hints pipeline, there are three steps described below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)
3. Configure the pipeline
    * [Local installation](configuration/local.md)
    * [Adding your own system](configuration/adding_your_own.md)

## 1) Install Miniconda2

In the RZcluster: 

`module load miniconda2` 

In the assembly cluster: 

`module load Miniconda2` 

Otherwise, install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/)

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `NF-hints` is specified as the pipeline name.

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/NF-hints/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/NF-hints-master
```

To stop nextflow from looking for updates online, you can tell it to run in offline mode by specifying the following environment variable in your ~/.bashrc file:

```bash
export NXF_OFFLINE='TRUE'
``` 