# Genome Annotation - Troubleshooting

## (RNA-seq) Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4.  If you are running Single end data make sure to specify `--singleEnd`

If the pipeline can't find your files then you will get the following error

```
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` wont work give you what you want Whilst `*{R1,R2}*.gz` will.

## RNA-seq consists of mixed library types - what to do?!

We found that, for novel genome projects, the generation of dedicated mRNA-seq data is generally recommended and a fairly typical step in the planning. However, 
you may be in a situation where you have to rely on "mixed" RNAseq data - say, stranded and unstranded reads, or paired versus single-end reads. Due to how ESGA is designed,
this is not easily supported. We want to make sure that the execution of the workflow is "mostly simple", but in order to accomodate a range of RNA-seq types, we'd have to force you to
use some sort of sample sheet to let the pipeline know about the various characteristics of your input data. For the time being, we have decided to not do that - but would consider adding something
like this later if there is demand for it. For now, this simply means that you won't be able to feed your "mixed" RNA-seq data into ESGA - but you could devise an assembly strategy with Trinity outside 
of our pipeline and feed the resulting FASTA file as EST (--ESTs) evidence. However, please note that Trinity too is not exactly designed to assembly mixed data like that. 

## Pipeline is very slow

The performance of the pipeline can be tuned in a number of ways. We have observed long run times for highly fragmented genomes or if the number of sequences
in blast and exonerate jobs is chosen too large. 

### Genome assembly 
Your genome assembly should, at most, contain thousands of scaffolds - ideally much less than that. Modern sequencing approaches such as linked-read sequencing
or long reads make this a feasible goal for many organisms. It is also advisable to exclude very short scaffolds (<5kb) from your assembly prior to annotation;
usually these will not contain any useful information but can increase the run time dramatically. 

### Size of evidence data
Another critical factor is of course the amount of sequence data that is used for annotation. Especially
the processing of RNA-seq reads can be very demanding on the compute infrastructure. So as a general recommendation, we suggest to work with "tens of
thousands of proteins", "hundreds of thousands of EST/Transcripts" and "tens of millions of RNA-seq reads".

## Deleting the pipeline work directory takes a very long time

This is a "known" issue with very deep directory trees under linux. The typical command `rm -Rf work` can then take hours, or days to complete. 

As an alternative, try this:

```bash
mkdir -p empty_dir
rsync -a --delete empty_dir/ work/
rm -Rf empty_dir work
```

## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us and open an issue here on github.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).


