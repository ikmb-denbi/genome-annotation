#!/usr/bin/env nextflow
/*
========================================================================================
                    IKMB - de.NBI | Genome Annotation Pipeline
========================================================================================
 Genome Annotation Pipeline. Started 2018-10-17.
 #### Homepage / Documentation
 https://git.ikmb.uni-kiel.de/m.torres/NF-hints.git
 #### Authors
 MTorres m.torres <m.torres@ikmb.uni-kiel.de> - https://git.ikmb.uni-kiel.de/m.torres>
 MHoeppner m.hoeppner <m.hoeppner@ikmb.uni-kiel.de> 
----------------------------------------------------------------------------------------
*/

// Make sure the Nextflow version is current enough
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nextflow_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please use a more recent version of Nextflow!\n" +
              "============================================================"
}

def helpMessage() {
  log.info"""
  =================================================================
   IKMB - de.NBI | Genome Annotation Pipeline | v${params.version}
  =================================================================
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run main.nf --genome 'Genome.fasta' --prots 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' -c config/slurm.config --nthreads 3

  Mandatory arguments:
  --genome		Genome reference
      
  At least one of:
  --proteins		Proteins from other species
  --ESTs		ESTs or transcriptome
  --reads		Path to RNA-seq data (must be surrounded with quotes)

  Options:
    -profile            Hardware config to use (optional, will default to 'standard')
    Programs to run:
    --trinity		Run transcriptome assembly with Trinity and produce hints from the transcripts [ true (default) | false ]
    --gth		Run GenomeThreader to produce hints from protein file [ true (default) | false ]
    --augustus		Run Augustus to predict genes [ true (default) | false ]
    --funAnnot		Run functional annotation using Annie [ true (default) | false ]
 	
    Programs parameters:
    --species		Species database for RepeatMasker [ default = 'mammal' ]
    --rm_lib		Additional repeatmasker library in FASTA format [ default = 'false' ]
    --model		Species model for Augustus [ default = 'human' ]
    --UTR		Allow Augustus to predict UTRs (results are not optimal and takes much longer) [ 'on' | 'off' (default) ]
    --iso		Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer) [ 'true' | 'false' (default) ]
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --uniprot		Fasta file with Uniprot proteins for functional annotation [ default = '/bin/Eumetazoa_UniProt_reviewed_evidence.fa' ]
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide Blast jobs [ default = 500 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nrepeats		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]
    --ninterpro		Chunks (# of sequences) to divide InterPro jobs [ default = 200 ]

    Other options:
    --singleEnd		Specifies that the input is single end reads [ true | false (default) ]
    --rnaseq_stranded	Whether the RNAseq reads were sequenced using a strand-specific method (dUTP) [ true | false (default) ]
    --outdir		The output directory where the results will be saved [ default = 'output' ]
    -name		Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
	helpMessage()
	exit 0
}

// -----------------------------
// Setting standard input values
// -----------------------------

Genome = file(params.genome)
if( !Genome.exists() ) exit 1; "No genome assembly found, please specify with --genome"

if (params.proteins) {
	Proteins = file(params.proteins)
	if( !Proteins.exists() ) exit 1, "Protein file not found: ${Proteins}. Specify with --proteins."
	println "Will run Exonerate and GenomeThreader on Protein file"
}

if ( params.ESTs ){
	ESTs = file(params.ESTs)
	if( !ESTs.exists() ) exit 1, "ESTs file not found: ${ESTs}. Specify with --ESTs."
	println "Will run Exonerate on EST/Transcriptome file"
}

if (params.reads){
	println "Found reads and will run Hisat2 on RNA-seq data"
}

if (params.rm_lib) {
	RM_LIB = file(params.rm_lib)
	if (!RM_LIB.exists() ) exit 1, "Repeatmask library does not exist (--rm_lib)!"
}

// Make it fail if basic requirements are unmet
if (!binding.variables.containsKey("Proteins") && !binding.variables.containsKey("ESTs") && params.reads == false) {
	exit 1, "At least one type of input data must be specified (--proteins, --ESTs, --reads)"
}

if (params.trinity == true && params.reads == false) {
	exit 1, "Cannot run Trinity de-novo assembly without RNA-seq reads (specify both --reads and --trinity)"
}

// give this run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

// ----------------------
// ----------------------
// Starting the pipeline
// ----------------------
// ----------------------
// Logic:
// - repeatmask the genome
// -- make blast database
// --- align proteins with blast
// ---- align proteins with exonerate
// --- align proteins with genomethreader
// --- align ESTs with Blast
// -- Align RNAseq reads with HiSat2
// -- Assembly Trinity transcripts (genome-guided)
// -- Generate hints from all available sources
// -- Run AUGUSTUS prediction

// --------------------------
// Set various input channels
// --------------------------

Channel.fromPath(Genome)
	.set { GenomeHisat }

// Split the genome for parallel processing
Channel
	.fromPath(Genome)
	.splitFasta(by: params.nrepeats, file: true)
	.set { FastaRM }

// if proteins are provided
if (params.proteins) {
        Channel
        .fromPath(Proteins)
        .splitFasta(by: params.nblast, file: true)
        .set { fasta_prots }
} else {
        fasta_prots = Channel.from(false)
        prot_exonerate_hints = Channel.create()
}

// if ESTs are provided
if (params.ESTs) {
        Channel
                .fromPath(ESTs)
                .splitFasta(by: params.nblast, file: true)
                .set {fasta_ests}
} else {
        fasta_ests = Channel.from(false)
        est_exonerate_hints = Channel.create()
}

// if GenomeThreader should be run
if (params.gth == false ) {
        fasta_prots_gth = Channel.from(false)
        gth_hints = Channel.create()
} else {
        Channel
        .fromPath(Proteins)
        .splitFasta(by: params.nblast, file: true)
        .set {fasta_prots_gth}
}

// ---------------------------
// RUN REPEATMASKER
//----------------------------

// generate a soft-masked sequence for each assembly chunk
process runRepeatMasker {

	tag "${chunk_name}"
	publishDir "${OUTDIR}/repeatmasker/chunks"

	input: 
	file(genome_fa) from FastaRM

	output:
	file(genome_rm) into FastaChunks

	script:
	// Provide a custom repeatmask database
	chunk_name = genome_fa.getName()

	options = ""
	if (params.rm_lib) {
		options = "-lib ${RM_LIB}"
	}
	genome_rm = genome_fa.getBaseName() + ".fa.masked"
	genome_gff = genome_fa.getBaseName() + ".fa.out.gff3"
	
	"""
		RepeatMasker -species ${params.species} -gff -xsmall  $options -q -pa ${task.cpus} $genome_fa	
	"""
}

// Merge the repeat-masked assembly chunks
process runMergeRMGenome {

        publishDir "${OUTDIR}/repeatmasker"

	input:
	file(fasta_chunks) from FastaChunks.collect()

	output:
	file(masked_genome) into RMtoBlastDB

	script:
	
	masked_genome = "${GENOME.getBaseName()}.rm.fa"
	"""
		cat ${fasta_chunks} > merged.fa
		fastasort -f merged.fa > masked_genome
		rm merged.fa
	"""	
}

// Turn genome into a masked blast database
// Generates a dust mask from softmasked genome sequence
process runMakeBlastDB {
	
	input:
	file(genome_fa) from RMtoBlastDB

	output:
	set file(db_nhr),file(db_nin),file(db_nsq) into blast_db_prots, blast_db_ests, blast_db_trinity

	script:
	dbName = genome_fa.baseName
	db_nhr = dbName + ".nhr"
	db_nin = dbName + ".nin"
	db_nsq = dbName + ".nsq"
	db_mask = dbName + ".asnb"
	target = file(db_nhr)
	
	"""
		convert2blastmask -in $genome_fa -parse_seqids -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out $db_mask
		makeblastdb -in $genome_fa –input_type blastdb -parse_seqids -mask_data $db_mask -dbtype nucl -out $dbName 
	"""
}

// ---------------------
// PROTEIN DATA PROCESSING
// ---------------------

// ----------------------------
// Protein BLAST against genome
// ----------------------------

process runBlastProteins {

	tag "${chunk_name}"
	publishDir "${OUTDIR}/evidence/proteins/tblastn/chunks"

	input:
	file(protein_chunk) from fasta_prots
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_prots

	output:
	file(protein_blast_report) into ProteinBlastReport

	script:
	db_name = blastdb_nhr.baseName
	chunk_name = protein_chunk.baseName
	protein_blast_report = "${chunk_name}.blast"

	"""
		tblastn -db $db_name -query $protein_chunk -db_soft_mask 30 -max_target_seqs 1 -outfmt 6 -num_threads ${task.cpus} > protein_blast_report
	"""

}

// Parse Protein Blast output

process Blast2QueryTargetProts {

	tag "ALL"
        publishDir "${OUTDIR}/evidence/proteins/tblastn/chunk/"

	input:
	file(blast_reports) from ProteinBlastReport.collect()

	output:
	file(query2target_result_uniq_prots) into query2target_uniq_result_prots
	
	script:
	query_tag = Proteins.baseName
	query2target_result_uniq_prots = "${query_tag}.targets"
	
	"""
	cat $blast_reports > merged.txt
	BlastOutput2QueryTarget.pl merged.txt 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq_prots
	"""
}
// split Blast hits for parallel processing

query2target_uniq_result_prots
	.splitText(by: params.nexonerate, file: true)
	.set{ query2target_chunk_prots }

// Run Exonerate on the blast regions

process runExonerateProts {

	tag "${query_tag}"
	publishDir "${OUTDIR}/evidence/proteins/exonerate/chunks", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk_prots
	
	output:
	file(exonerate_chunk) into exonerate_result_prots
	
	script:
	query_tag = Proteins.baseName
	exonerate_chunk = "${hits_chunk.baseName}.${query_tag}.exonerate.out"
	
	"""
		runExonerate_fromBlastHits_prot2genome.pl $hits_chunk $Proteins $Genome
		grep -v '#' exonerate.out | grep 'exonerate:protein2genome:local' > $exonerate_chunk
	"""
}

process runMergeExonerateHints {

	publishDir "${OUTDIR}/evidence/proteins/exonerate/"

	input:
	file(chunks) from exonerate_result_prots.collect()

	output:
	file(exonerate_hints) into prot_exonerate_hints

	script:
	exonerate_gff = "proteins.exonerate.hints.gff"
	"""
		cat $chunks > all_chunks.out
		Exonerate2GFF_protein.pl all_chunks.out $exonerate_gff
	"""
}


// ------------------------------------
// GenomeThreader hints generation
// ------------------------------------

// Run genome threader for proteins if requested
process runGenomeThreaderProteins {

	tag "${query_tag}"
	publishDir "${OUTDIR}/evidence/proteins/gth/chunks/"

	input:
	file(protein_chunk) from fasta_prots_gth

	output:
	file(gth_chunk) into ProteinGTHChunk

	script:
	query_tag = protein_chunk.getName()
	protein_tag = Proteins.baseName
	gth_chunk = "${query_tag}.${protein_tag}.gth"

	"""
		gth -genomic $Genome -protein $protein_chunk -gff3out -intermediate -o $gth_chunk
	"""	
}

process GenomeThreader2HintsProts {

	tag "${query_tag}"
        publishDir "${OUTDIR}/evidence/proteins/gth/chunks"

	input:
	file(not_clean_gth) from ProteinGTHChunk
	
	output:
	file(gth_hints) into ProteinGTHChunkHint
	
	script:
	gth_hints = not_clean_gth.baseName + ".clean"
	
	"""
	gt gff3 -addintrons yes -setsource gth -tidy yes -addids no $not_clean_gth > not_clean_gth_wIntrons
	grep -v '#' not_clean_gth_wIntrons > no_hash_gth
	GTH_rename_splitoutput.pl no_hash_gth > clean_gth
	grep -e 'CDS' -e 'exon' -e 'intron' clean_gth | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > gth_hints
	
	"""
}

process GenomeThreaderMergeHints {

        publishDir "${OUTDIR}/evidence/proteins/gth/"

	input:
	file(gth_hint_chunks) from ProteinGTHChunkHint.collect()

	output:
	file(merged_gth_hints) into gth_hints

	script:
	gth_hints = " gth.proteins.hints.gff"

	"""
		cat $gth_hint_chunks >> merged_gth_hints
	"""
}

// --------------------------
// -------------------
// EST DATA PROCESSING
//-------------------
// --------------------------


/*
 * EST blasting
*/

// Blast the ESTs against the nucleotide database

process runBlastEst {

	publishDir "${OUTDIR}/evidence/EST/blast/chunks"

	input:
	file(est_chunk) from fasta_ests
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_ests
		
	output:
	file(blast_report) into ESTBlastReport

	when:
	params.ESTs

	script:
	db_name = blastdb_nhr.baseName
	chunk_name = est_chunk.baseName

	"""
		blastn -db $db_name -query $query_fa_ests -max_target_seqs 1 -outfmt 6 -num_threads ${task.cpus} > blast_report
	"""
}

// Parse the EST Blast output
process Blast2QueryTargetEST {

        publishDir "${OUTDIR}/evidence/EST/blast"

	input:
	file(blast_report) from ESTBlastReport.collectFile()

	output:
	file(target_list) into query2target_uniq_result_ests

	script:
	query_tag = ESTs.baseName
	query2target_result_uniq_ests = "EST.blast.targets.txt"

	"""
	BlastOutput2QueryTarget.pl $blast_report 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq_ests
	"""

}

query2target_uniq_result_ests
	.splitText(by: params.nexonerate, file: true).set{query2target_chunk_ests}	

// Run exonerate on the EST Blast chunks
process runExonerateEST {

	tag "${query_tag}"
	publishDir "${OUTDIR}/evidence/EST/exonerate/chunks", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk_ests
	
	output:
	file(results) into exonerate_result_ests
	
	script:
	query_tag = ESTs.baseName
	chunk_name = hits_chunk.getName()
	results = "${chunk_name}.${query_tag}.exonerate.out"
	
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_chunk $ESTs $Genome
	mv exonerate.out $results
	"""
}

// Convert exonerate hits to hints

process Exonerate2HintsEST {

	tag "${query_tag}"
	publishDir "${OUTDIR}/evidence/EST/exonerate/", mode: 'copy'

	input:
	file(exonerate_result) from exonerate_result_ests.collectFile()
	
	output:
	file(exonerate_gff) into output_gff_ests
	
	script:
	query_tag = ESTs.baseName
	exonerate_hints = "${query_tag}.EST.hints.gff"
		
	"""
	grep -v '#' $exonerate_result_ests | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_EST.pl exonerate_gff_lines exonerate_gff
	"""
}
