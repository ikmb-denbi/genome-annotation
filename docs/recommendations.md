# Recommendations

Evidence data for annotation can be obtained from a number of sources. Below follow some general suggestions that we found to work well.

## Proteins

[Uniprot](https://www.uniprot.org/) is generally a good place to start. The search interface allows you to limit results by taxnomic group and level of experimental support for 
each sequence. As a rule of thumb, try to focus on sequences that are flagged as "full length" and that have experimental support from either
protein or transcriptome sequencing. 

[RefSeq](https://www.ncbi.nlm.nih.gov/protein/) is another useful database. However, sequences here tend to include a lot more computational predictions,
which runs the risk of perpetuating annotation errors. 

## ESTs

EST data can be downloaded from [GenBank](https://www.ncbi.nlm.nih.gov/nucleotide) or [ENA](https://www.ebi.ac.uk/ena). 
Make sure that the sequences are from your species of interest; on the nucleotide level even relatively small evolutionary timescales could
greatly dimminish the alignment quality/rate. 

An alternative to traditional ESTs are sequence data from de-novo assembled transcriptomes. For this, we recommend genome-guided [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly).
Please be advised that this pipeline actually performs this step for you, if you specify both `--reads` and `-trinity true`.

## RNA-seq data

This is perhaps the most important type of annotation evidence, as it stems specifically from your organism of interest and should provide
a good resolution of splice junctions. 

We recommend poly-A selected, strand-specific paired-end RNA-seq libraries for the purpose of annotation (any length is fine, but longer reads
tend to give better results for de-novo transcriptome assembly). Ideally, these libraries cover different develpomental stages and body parts.
For mammals, a typical sequencing depth per library is around 25 million PE reads (this is more related to the size of the transcriptome than 
of the genome!)

## Repeats

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
