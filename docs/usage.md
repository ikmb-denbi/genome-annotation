# How to run the pipeline 

![](../images/genome-annotation_dag_mod.svg) 

The typical command for running the pipeline is as follows:

`nextflow run main.nf --genome 'Genome.fasta' --prots 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' --ESTs 'ESTs.fa' --nthreads 3 --outdir 'my_species_annotation_out` 
  
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

#### `--species` [ default = 'mammal' ]
Species database for RepeatMasker. 

#### `--model` [ default = 'human' ]
Species model for Augustus. 

#### `--UTR` [ 'on' | 'off' (default) ] 
Allow Augustus to predict UTRs (results are not optimal and takes much longer). 

#### `--isof` [ 'true' | 'false' (default) ] 
Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer). 

#### `--augCfg` [ default = 'bin/augustus_default.cfg' ]
Location of augustus configuration file. 

#### `--uniprot` [ default = '/bin/Eumetazoa_UniProt_reviewed_evidence.fa' ]
Fasta file with Uniprot proteins for functional annotation. 
    
### 4. How to split programs 
#### `--nblast` [ default = 500 ]
Chunks (# of sequences) to divide Blast jobs. 

#### `--nexonerate` [ default = 200 ]
Chunks (# of blast hits) to divide Exonerate jobs.

#### `--nrepeats` [ default = 30 ]
Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs.

#### `--ninterpro` [ default = 200 ] 
Chunks (# of sequences) to divide InterPro jobs. 

#### `--nthreads` [default = 1] 
Number of cpus for programs that allow multi-threaded mode. 

### 5. Other options 
#### `--singleEnd` [ true | false (default) ]
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run. 

#### `--outdir` [ default = 'Hints_annotation_output' ]
The output directory where the results will be saved. 

#### `--allHints` [ default = 'AllHints.gff' ] 
Name of final GFF file with all hints. 

#### `--addHints` [ default = 'false' ]
Additional hints file (in GFF format), to be concatenated to the resulting hints before running augustus. 

#### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different combination of compute environment and installation estrategy (see [Installation instructions](../docs/Installation.md)). 

Available profiles are:

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
