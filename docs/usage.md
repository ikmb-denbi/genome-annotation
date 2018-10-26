# How to run the pipeline 

![](../images/genome-annotation_dag_mod.svg) 

### 1. Mandatory arguments 

#### `--genome` 
Location of the genome you want to annotate. It must be in FASTA format. 

### 2. Evidences. At least one of:

#### `--reads` 
Location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

#### `--ESTs` 
Location of a single FASTA file with all EST sequences or assembled transcriptome(s) from the species of interest. If you have multiple files, concatenate them before into a single file. 

#### `--prots` 
Location of a single FASTA file with protein sequences from related species. If you have multiple files, concatenate them before into a single file. 

### 3. Programs to run 
By default, the complete pipeline you see above will run. You can skip some steps if you already have created hints files (see `--addHints`), if you already have assembled a transcriptome, if you don't want to run gene prediction and/or functional annotation, if you don't want to install some of the required programs or if you want the pipeline to finish faster. 

#### `--trinity` [ true (default) | false ] 
Run transcriptome assembly with Trinity and produce hints from the transcripts. 

#### `--gth` [ true (default) | false ] 
Run GenomeThreader to produce hints from protein file.  

#### `--RM` [ true (default) | false ] 
Run RepeatMasker to produce hints. 

#### `--augustus` [ true (default) | false ] 
Run Augustus to predict genes.  

#### `--funAnnot` [ true (default) | false ] 
Run functional annotation using Annie. 

### 5. Parameters for specific programs  

Programs parameters:
    --species		Species database for RepeatMasker [ default = 'mammal' ]
    --model			Species model for Augustus [ default = 'human' ]
    --UTR			Allow Augustus to predict UTRs (results are not optimal and takes much longer) [ 'on' | 'off' (default) ]
    --isof			Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer) [ 'true' | 'false' (default) ]
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --uniprot		Fasta file with Uniprot proteins for functional annotation [ default = '/bin/Eumetazoa_UniProt_reviewed_evidence.fa' ]
    
### 4. How to split programs 


    
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide Blast jobs [ default = 500 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nrepeats		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]
    --ninterpro		Chunks (# of sequences) to divide InterPro jobs [ default = 200 ]
    --nthreads		Number of cpus for programs that allow multi-threaded mode [default = 1]	

### 5. Other options 
    Other options:
    --singleEnd		Specifies that the input is single end reads [ true | false (default) ]
    --outdir		The output directory where the results will be saved [ default = 'Hints_annotation_output' ]
    --allHints		Name of final GFF file with all hints [ default = 'AllHints.gff' ]
    --addHints		Additional hints file (in GFF format), to be concatenated to the resulting hints before running augustus [ default = 'false' ]
    -name			Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

#### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

#### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment - follow the links below to see instructions for running on that system. Available profiles are:

* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Runs using the `local` executor and pulls software from dockerhub: [`NF-hints`](http://hub.docker.com/r/NF-hints/)
* `aws`
    * A starter configuration for running the pipeline on Amazon Web Services. Uses docker and Spark.
    * See [`docs/configuration/aws.md`](configuration/aws.md)
* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### 6. Nextflow parameters

#### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

#### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

#### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

#### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

#### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1` 

### Job Resources
#### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

#### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples. 

## All the commands at a glance 
```
The typical command for running the pipeline is as follows:

  nextflow run main.nf --genome 'Genome.fasta' --prots 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' -c config/slurm.config --nthreads 3

  Mandatory arguments:
  --genome		Genome reference
  -profile		Hardware config to use
      
  At least one of:
  --prots		Proteins from other species
  --ESTs		ESTs or transcriptome
  --reads		Path to RNA-seq data (must be surrounded with quotes)

  Options:
    Programs to run:
    --trinity		Run transcriptome assembly with Trinity and produce hints from the transcripts [ true (default) | false ]
    --gth			Run GenomeThreader to produce hints from protein file [ true (default) | false ]
    --RM			Run RepeatMasker to produce hints [ true (default) | false ]
    --augustus		Run Augustus to predict genes [ true (default) | false ]
    --funAnnot		Run functional annotation using Annie [ true (default) | false ]
 	
    Programs parameters:
    --species		Species database for RepeatMasker [ default = 'mammal' ]
    --model			Species model for Augustus [ default = 'human' ]
    --UTR			Allow Augustus to predict UTRs (results are not optimal and takes much longer) [ 'on' | 'off' (default) ]
    --isof			Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer) [ 'true' | 'false' (default) ]
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --uniprot		Fasta file with Uniprot proteins for functional annotation [ default = '/bin/Eumetazoa_UniProt_reviewed_evidence.fa' ]
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide Blast jobs [ default = 500 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nrepeats		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]
    --ninterpro		Chunks (# of sequences) to divide InterPro jobs [ default = 200 ]
    --nthreads		Number of cpus for programs that allow multi-threaded mode [default = 1]	

    Other options:
    --singleEnd		Specifies that the input is single end reads [ true | false (default) ]
    --outdir		The output directory where the results will be saved [ default = 'Hints_annotation_output' ]
    --allHints		Name of final GFF file with all hints [ default = 'AllHints.gff' ]
    --addHints		Additional hints file (in GFF format), to be concatenated to the resulting hints before running augustus [ default = 'false' ]
    -name			Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
```
