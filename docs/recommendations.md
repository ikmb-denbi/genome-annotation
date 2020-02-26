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
Please be advised that this pipeline actually performs this step for you, if you specify both `--reads` and `--trinity true`.

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

## Training a new Augustus model

Augustus ships with a number of high-quality prediction models from a range of taxonomic groups. Usually, the easiest approach is thus to
use a model that is taxonomically somewhat close to the species you are trying to annotate. For example, if annotating a bird, the built-in model for
chicken should work just fine. Remember, Augustus uses a range of hints produced by this pipeline to inform its gene finding anyway. As long as the basics of
the model are appropriate for you organism, this is a good approach.

However, if you are not getting satisfying results or find that no pre-existing model is likely appropriate for your genome of choice, you can enable an
automatic training routine. This will use available transcriptome data to either refine an existing model (--model exists) or built one from scratch 
(--model does not yet exist). This is somewhat experimental and depends a lot on the quality of the input data. It also takes a pretty long time (several
days for larger genomes) as it needs to first re-construct gene models from the aligned transcriptome data, select all the models that are probably
full length and finally use these gene structures to train Augustus. 

## Improving an initial gene build

Chances are the first run of the pipeline will not produce satisfactory results. Apart from trying to add additional data (e.g. more RNAseq, or ESTs, etc), an obvious
aspect is the tuning of parameters for AUGUSTUS and EvidenceModeler. Both tools require files that control their behavior. 

Changes to these files will force parts of the pipeline to rerun that rely on stages that make use of these files. EVM and AUGUSTUS are at the very end of the ESGA workflow, so resumeing the pipeline with
modified config files typically does not take too long. You can thus try to improve your annotation iteratively by just playing with these files. More below. 

### Augustus
AUGUSTUS uses extrinsic hints to guide the gene finding. These hints are produced by the ESGA pipeline from provided evidence. How these hints are weighted inside AUGUSTUS however are 
controlled by a config file. 

Please see [here](https://github.com/Gaius-Augustus/Augustus/blob/master/config/extrinsic/extrinsic.cfg) for instructions on how to tweak these parameters for optiomal performance. ESGA uses, by default,
an extrinsic config file that we have set up to work for our typical projects. You can however pass a modified version from the command line instead using the `--augCfg` option together with `-resume`. 

## EVM
ESGA produces several types of inputs (all of them optional) that are then combined by EVM into a consensus gene buid. The weight given to each type of input is controlled by the [weights.txt](../assets/evm/weights.txt) file. 
EVM uses weights between 1 and 10 to determine which evidence to consider first. To pass your own, modified version of the weights file, use the `--evm_weights my_weights.txt` syntax together with the `-resume` option. 
