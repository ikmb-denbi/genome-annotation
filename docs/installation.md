# Installation and configuration 

## Installing Nextflow 

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
Therefore the first thing to do is install Nextflow. 

If you are working in the **IKMB RZ cluster**, you can simply load the following modules:

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

**You need NextFlow version >= 18.01.0 to run this pipeline.** 

## Cloning the genome-annotation repository 

To run the pipeline you first have to check out the code to a location in your system (i.e. $HOME/git/). 

``` 
cd $HOME/git/ 

git clone git@github.com:ikmb-denbi/genome-annotation.git
``` 
 
## Installing all other software 

This pipeline uses a lot of different bioinformatics software (you can see a full list at the end of this document). How you proceed to install these programs will depend on what system/cluster you are using: 

### A. Working in the IKMB RZ cluster 

In the **IKMB RZ cluster**, all these programs are available as modules and will be loaded automatically as the pipeline runs. You don't need to do anything else, just make sure that you run the pipeline using the parameter `-profile standard`. This is anyway the default profile, so you don't even need to especify it, only don't use any other. 

### B. Not in the IKMB? Use Bioconda 

Most of the required programs are available as [bioconda packages](https://bioconda.github.io/recipes.html) for easy installation. All you need to do is install the corresponding miniconda2 for your system: 

[miniconda2 installer](https://repo.continuum.io/miniconda/) 

In this case, when you run the pipeline add `-profile conda` (if in the IKMB) or `-profile conda_custom` (see below) to your command. 

Nextflow will install the environments with the necessary packages as it runs. 

#### Missing programs: 

Some of the required programs are not available as conda packages yet, so you will have to install them yourself. Some parts of the pipeline will run successfully anyway, but you need to turn off the ones that don't. 

1. **GenomeThreader:** it is used to map the protein evidences to the genome and create hints. It is not necessary to run the pipeline and, by default, protein evidence will be mapped using [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate). However, my experience shows that the two programs complement each other: some models that were missed by Exonerate are found by GenomeThreader, and viceversa. You can download GenomeThreader from [here](http://genomethreader.org/download.html) and follow the installation instructions. If you want to run the pipeline without this program, you must use `--gth false` on your command call.

2. **Annie:** it is used to transfer the functional annotations to your gff3 file. You can download it from [here](http://genomeannotation.github.io/annie/) and extract it in your system. Make sure the executable file `annie.py` is in your path. If you don't want to perform functional annotation, use `--funAnnot false`. 

3. **InterProScan:** it is also necessary during the functional annotation step. It it is used to scan the predicted protein sequences for known protein signatures (functional domains, GO terms, etc) searching in different public databases. You can download it from [here](https://www.ebi.ac.uk/interpro/download.html) and extract it in your system. Make sure the executable `interproscan.sh` is in your path. If you want to run the pipeline without this program, use `--funAnnot false`. 

4. **Bioruby:** to run the functional annotation step you also need to have the bioruby library installed. You can download and install ruby if you haven't done so yet from [here](https://www.ruby-lang.org/en/). Then install bioruby using the RubyGems tool: `gem install bio` 

*coming soon:* [Singularity]() image with these programs installed. 

#### Configuration for other clusters (`conda_custom`)

It is entirely possible to run this pipeline on clusters other than the IKMB, but you will need to set up your own config file so that the pipeline knows how to work with your cluster. 

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline. 

If you are the only person to be running this pipeline, you can simply modify the [custom.config](../conf/custom.config) file to include the especifications of your cluster. Then remember to run the pipeline with `-profile conda_custom`. This profile automatically uses the [custom.config](../conf/custom.config) file, so keep the file name! 

If you are familiar with Nextflow, you can add your own profiles by changing the [nextflow.config](../nextflow.config) file.

#### Configuration for local execution (`local`)

If you don't want to use any cluster management system but rather run the pipeline locally, use `-profile local`. 

### C. Install all programs yourself 

Here is a list of all the programs necessary to run the complete genome-annotation pipeline (`--prots proteins.fa --ESTs ESTs.fa --reads '*_R{1,2}.fastq' --gth true --RM true --trinity true --augustus true --funAnnot true`).

1. [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.2.30 

2. [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) v2.2.0 

3. [Bioperl](https://bioperl.org) 

4. [GenomeThreader](http://genomethreader.org/download.html) v1.7.0 

5. [Genometools](http://genometools.org) v1.5.6 

6. [RepeatMasker](http://www.repeatmasker.org) v.4.0.6 

7. [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.4.4 

8. [Hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml) v2.1.0 

9. [Samtools](http://www.htslib.org/download/) v1.5 

10. [Augustus](http://bioinf.uni-greifswald.de/augustus/) v3.2.1 

11. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) v2.5.1 

12. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.2.3 

13. [Ruby](https://www.ruby-lang.org/en/) v2.2.2 + bioruby library 

14. [Annie](http://genomeannotation.github.io/annie/)  

15. [Interproscan](https://www.ebi.ac.uk/interpro/download.html)  

You can use the `-profile self_install` to run the pipeline in this case. Modify the **_genome-annotation/conf/custom.config_** as necessary for your system. The pipeline has been tested successfully with the versions that are here described. 

