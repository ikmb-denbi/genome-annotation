# How to run the pipeline 

The typical command for running the pipeline is as follows:

```
nextflow run main.nf --genome 'Genome.fasta' --proteins 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' --ESTs 'ESTs.fa' --outdir 'my_species_annotation_out'
```

This will run all the steps in the pipeline (Proteins, ESTs/Transcriptome, RNAseq). The types of evidences you provide determine which parts of the pipeline are actually run. 


### Parameters file

In the next section, you will find a list of all user-configurable pipeline options. 
You can of course provide each option as a command line parameter. But this can get a bit tedious. As an alterantive, you can provide a configuration file using the YAML format. An example is included under [../assets/config.yaml](../assist/config.yaml). To provide a config file as an option, 
use `-params-file my_config.yaml`. The revised command could then read:

`nextflow run /path/to/main.nf -params-file config.yaml -profile your_profile`

The full YAML options file:

```yaml
genome: "/path/to/genome.fa"
proteins: "/path/to/proteins.fa"
ESTs: "/path/to/ests.fa"
rm_lib: "/path/to/repeats.fa"
rm_species: "mammal"
reads: "/path/to/*_R{1,2}_001.fastq.gz"
trinity: false
augustus: true
model: "human"
training: false
augCfg: false
nblast: 200
blast_evalue: 0.001
nexonerate: 100
nrepeats: 20
singleEnd: false
rnaseq_stranded: false
outdir: "output"
run_name: "evidences"
max_intron_size: 200000
```
  
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
Location of a single FASTA file with all EST sequences or assembled transcriptome(s) from the species of interest. If you have multiple files, concatenate 
them into a single file first and make sure that the sequence names are not duplicated (this happens when you try to merge two Trinity assemblies, 
for example). 

#### `--proteins` 
Location of a single FASTA file with protein sequences from related species. If you have multiple files, concatenate them before into a single file. 

### 3. Programs to run 
By default, the complete pipeline you see above will run, given the types of evidences your provide. However, you can skip some steps if you want. 
For example,if you already have assembled a transcriptome or if you don't want to run gene prediction. 

#### `--trinity` [ true (default) | false ] 
Run transcriptome assembly with Trinity and produce hints from the transcripts. 

#### `--augustus` [ true (default) | false ] 
Run Augustus to predict genes.  

### 4. Within-scaffold parallelization

#### `--max_intron_size <int>` [ 20000 (default) ]
This pipeline will run certain programs not on full scaffolds, but on clusters of data within those scaffolds. The factor determining how to build these evidence clusters is the expected maximum intron size for your organism of interest. 
The default value is set to 20000 - for something like a nematode, this would be too long, for human it would probably be fine, although 
a few introns are much longer. Genes containing such extraordinarily large introns will then probably be mis-annotated. 
However, choosing too large values will drastically increase the run time. 

### 5. Parameters for specific programs 
To run some of the programs, additional information is required. All options have a (hopfully reasonable) default, but you must check if it is the proper one for your organism and for the output you expect. 

#### `--rm_species` [ default = 'mammal' ]
Species database for RepeatMasker. This option will use the very limited built-in database of Repeatmasker and is only useful for annotating primates (mostly). It is also mutually exclusive with the preferred option `--rm_lib <file>` (below). 

#### `--rm_lib`[ fasta file | false ]
By default, Repeatmasker will run with the built-in DFam hmm profile for (mostly) primates. It is thus generally advisable to instead provide 
repeat annotations in FASTA format. Possible sources include self-computed repeats (using RepeatModeler) or curated repeat libraries from 
GRINST (www.grinst.org, commercial). 
If you have a copy of the complete Repeatmasker library (and an installation of RM), you can extract the repeat annotation from a species like this: 

```
# Get the Tree of all available species: 
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -tree

# Select the name of the species (e.g. "Ostreoida") from the output tree and do:
perl /[...]/RepeatMasker/4.0.8/util/queryRepeatDatabase.pl -species Ostreoida > RMdb_Ostreoida.fa
``` 

Then run the pipeline with the option "--rm_lib RMdb_Ostreoida.fa". 

Please note: If *no* repeats are specified, the pipeline will try to model repeats de-novo using the RepeatModeler package. However, this requires a
sufficient number of related repeats to be present in your assembly. If your genome was assembled from short reads, this strategy may not return anything -
short read assemblies tend to collapse repeats. In this case, the pipeline will fall back to the built-in library that ships with RepeatMasker
(which is very limited and probably not very useful).

#### `--model` [ default = 'human' ]
Species model for Augustus. A list of valid identifiers can be found [here](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md).

#### `--training` [ true | false (default) ]
Perform training of an Augustus prediction profile. Use this option if you want to fine-tune an existing model (--model specifies an existing Augustus profile)
or train a new model from scratch (--model does not specify an existing Augustus profile). This requires the RNA-seq gene builder PASA to be run, so you
have to provide either RNA-seq reads (and --trinity) or a reasonably sized set of ESTs/assembled transcripts (or both).

#### `--training_models [ default = 1000 ]
Defines how many of the gene models reconstructed by PASA and flagged as "complete" are to be used in the training of AUGUSTUS. The default of 1000 is usually fine. 

#### `--augCfg` [ default = 'bin/augustus_default.cfg' ]
Location of Augustus configuration file. By default, this pipeline uses config file that we found to work well for predicting gene models in mammalian genomes using the kinds of extrinsic hints constructed by this pipeline. 
    
### 6. How to tune the speed of the pipeline - data splitting

One of the advantages of using Nextflow is that it allows you to speed up a pipeline by splitting some of the input files into smaller chunks before 
running specific programs. Then that program can be run on each smaller chunk in parallel in a compute cluster. 
When all instances of the program are finished, Nextflow can correctly put together all the results in a single output for that program. Depending on the size and contiguity of your target genome and the size of the evidence data, you may want to tweak one or several of the parameters below. If unsure, 
leave at the defaults.

#### `--chunk_size` [ default = 50000 ]
Length of sequence to which break the genome down to for faster Blastx searches. Longer sequences will increase the runtime. This parameter is directly related to `--nblast` - which defines how many sequences are put into each blastx search.

#### `--nblast` [ default = 500 ]
Number of sequences in each Blast job. Larger values will usually create longer run times, but decrease the number of parallel jobs and load on the file system. 

#### `--nexonerate` [ default = 200 ]
Number of alignments to compute in each Exonerate job. Larger values will usually create longer run times, but decrease the number of parallel jobs and load on the file system.

#### `--nchunks` [ default = 10 ]
Number of pieces to split the genome assembly into before running RepeatMasker and Augustus jobs. For example, if you want to annotate the human genome, a chunk size of 23 would place each chromosome into one chunk. Chunks cannot be smaller than individual scaffolds/chromosomes.

### 7. Other options 

#### `--email` [ you@somewhere.com | false (default)]
If you specify an Email address, the pipeline will send a notification upon completion. However, for this to work, the node running the nextflow process must have a configured Email server. 

#### `--singleEnd` [ true | false (default) ]
By default, the pipeline expects paired-end RNA-seq data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run. 

#### `--outdir` [ default = 'annotation_output' ]
The output directory where the results will be saved. 

#### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different combination of compute environment and installation estrategy (see [Installation instructions](../docs/installation.md)). 

### 8. Nextflow parameters (indicate with single dash "-")

#### `-params-file config.yaml`
All the above options can be passed from either the command line or through a configuration file. A suitable template is included under assets/config.yaml.

#### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

#### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

### Job Resources
#### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

#### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples. 

