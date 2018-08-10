#!/usr/bin/env nextflow
/*
========================================================================================
                         NF-Hints
========================================================================================
 NF-hints Analysis Pipeline. Started 2018-08-03.
 #### Homepage / Documentation
 https://git.ikmb.uni-kiel.de/m.torres/NF-hints.git
 #### Authors
 MTorres m.torres <m.torres@ikmb.uni-kiel.de> - https://git.ikmb.uni-kiel.de/m.torres>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     NF-hints v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run Protein2Hints.nf --genome 'Genome.fasta' --query 'Proteins.fasta' -profile docker

    Mandatory arguments:
      --genome                      Genome reference
      -profile                      Hardware config to use. docker / aws
      
    At least one of:
      --prots						Proteins from other species
      --ESTs						ESTs or transcriptome
      --reads						Path to RNA-seq data (must be surrounded with quotes)

    Options:
	  --trinity						Run transcriptome assembly with Trinity and produce hints from the transcripts [ true (default) | false ]
	  --gth							Run GenomeThreader to produce hints from protein file [ true (default) | false ]
	  --RM							Run RepeatMasker to produce hints [ true (default) | false ]
      --nblast						Chunks to divide Blast jobs [ default = 10 ]
      --nexonerate					Chunks to divide Exonerate jobs [ default = 10 ]
	  --nrepeats					Chunks to divide RepeatMasker jobs [ default = 2 ]
	  --nthreads					Number of cpus for programs that allow multi-threaded mode [default = 1]
	  --species						Species database for RepeatMasker [ default = 'mammal' ]
	  --singleEnd                   Specifies that the input is single end reads [ true | false (default) ]

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}



//Default variables:
params.prots = false
params.ESTs = false
params.reads = false

params.trinity = true
params.gth = true
params.RM = true

params.nblast = 10
params.nexonerate = 2
params.nrepeats = 2
params.nthreads = 1
params.species = "mammal"
params.name = false
params.singleEnd = false


// Validate inputs
if ( params.genome ){
	Genome = file(params.genome)
    if( !Genome.exists() ) exit 1, "Genome file not found: ${Genome}"
}

x = 0

if ( params.prots ){
	Proteins = file(params.prots)
	x = x + 1
    if( !Proteins.exists() ) exit 1, "Protein file not found: ${Proteins}"
    println "Will run Exonerate and GenomeThreader on Protein file"
}

if ( params.ESTs ){
	ESTs = file(params.ESTs)
	x = x + 1
    if( !ESTs.exists() ) exit 1, "ESTs file not found: ${ESTs}"
    println "Will run Exonerate on EST/Transcriptome file"
}

if (params.reads){
	x = x + 1
	println "Will run Hisat2 on RNA-seq data"
}

if (x == 0) { 
	exit 1, "At least one data file must be especified"
}


if (params.trinity == true && params.reads == false) {
	exit 1, "Cannot run Trinity without RNA-seq reads"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Header log info
log.info """=======================================================
                                            
    ___  ___          __   __   __   ___     
    |__| |__  \\/  __ /  ` /  \\ |__) |__        
    |__| |    /\\     \\__, \\__/ |  \\ |___    
                                               

NF-hints v${params.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'NF-hints'
summary['Pipeline Version'] = params.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Fasta Ref']    = Genome
summary['Proteins']		= params.prots
summary['ESTs']			= params. ESTs
summary['Reads']		= params.reads
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}



Channel
	.fromPath(Genome)
	.set { inputMakeblastdb }
	
// We check if the blast db already exists - if not, we create it

/*
 * STEP 1 - Make Blast DB
 */
 
process RunMakeBlastDB {
	
	tag "${dbName}"
	publishDir "${params.outdir}/BlastDB", mode: 'copy'
	
	input:
	file(genome) from inputMakeblastdb
	
	output:
	set file(db_nhr),file(db_nin),file(db_nsq) into blast_db_prots, blast_db_ests, blast_db_trinity
	
	script:
	dbName = genome.baseName
	db_nhr = dbName + ".nhr"
	db_nin = dbName + ".nin"
	db_nsq = dbName + ".nsq"

	target = file(db_nhr)
	
    if (!target.exists()) {
		"""
		makeblastdb -in $genome -dbtype nucl -out $dbName
		"""
	}
	
}


/*
 * Proteins Block
 */

// Create a channel emitting the query fasta file(s), split it in chunks 

if (params.prots) {
	Channel
		.fromPath(Proteins)
		.splitFasta(by: params.nblast, file: true)
		.set {fasta_prots}
} else { 
	fasta_prots = Channel.from(false)
}


/*
 * STEP Proteins.1 - Blast
 */
 
process RunBlastProts {
	
	tag "${chunk_name}"
	publishDir "${params.outdir}/blast_results/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa_prots from fasta_prots 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_prots.collect()
	
	output:
	file blast_result_prots
	
	when:
	params.prots != false
		
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa_prots.baseName
	
	"""
	tblastn -db $db_name -query $query_fa_prots -max_target_seqs 1 -outfmt 6 -num_threads $params.nthreads > blast_result_prots
	"""
}


/*
 * STEP Proteins.2 - Parse Blast Output
 */

process Blast2QueryTargetProts {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/blast2targets", mode: 'copy'
	
	input:
	file all_blast_results_prots from blast_result_prots.collectFile()
	
	output:
	file query2target_result_uniq_prots into query2target_uniq_result_prots
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	"""
	BlastOutput2QueryTarget.pl $all_blast_results_prots 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq_prots
	"""
}

query2target_uniq_result_prots
	.splitText(by: params.nexonerate, file: true).set{query2target_chunk_prots}	


/*
 * STEP Proteins.3 - Exonerate
 */
 
process RunExonerateProts {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/exonerate/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk_prots
	
	output:
	file 'exonerate.out' into exonerate_result_prots
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	runExonerate_fromBlastHits_prot2genome.pl $hits_chunk $Proteins $Genome
	"""
}


/*
 * STEP Proteins.4 - Exonerate to Hints
 */
 
process Exonerate2HintsProts {
	
	tag "${query_tag}"
	
	input:
	file exonerate_result_prots
	
	output:
	file exonerate_gff into output_gff_prots
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	grep -v '#' $exonerate_result_prots | grep 'exonerate:protein2genome:local' > exonerate_gff_lines
	Exonerate2GFF_protein.pl exonerate_gff_lines exonerate_gff
	"""
}

output_gff_prots
 	.collectFile(name: "${params.outdir}/Hints/Hints_proteins_exonerate.gff")


/*
 * STEP Proteins.5 - GenomeThreader
 */
 
process RunGenomeThreaderProts {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/genomethreader", mode: 'copy'
		
	output:
	file output_gth
	
	when:
	params.prots != false && params.gth != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	gth -genomic $Genome -protein $Proteins -gff3out -intermediate -o output_gth
	"""
}


/*
 * STEP Proteins.6 - GenomeThreader to Hints
 */
 
process GenomeThreader2HintsProts {
	
	tag "${query_tag}"
	
	input:
	file not_clean_gth from output_gth
	
	output:
	file gth_hints
	
	when:
	params.prots != false && params.gth != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	gt gff3 -addintrons yes -setsource gth -tidy yes -addids no $not_clean_gth > not_clean_gth_wIntrons
	grep -v '#' not_clean_gth_wIntrons > no_hash_gth
	GTH_rename_splitoutput.pl no_hash_gth > clean_gth
	grep -e 'CDS' -e 'exon' -e 'intron' clean_gth | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > gth_hints
	"""
}

gth_hints
	.collectFile(name: "${params.outdir}/Hints/Hints_proteins_genomethreader.gff")



/*
 * ESTs Block
 */
 
if (params.ESTs) {
	Channel
		.fromPath(ESTs)
		.splitFasta(by: params.nblast, file: true)
		.set {fasta_ests}
} else { 
	fasta_ests = Channel.from(false)
}


/*
 * STEP ESTs.1 - Blast
 */
 
process RunBlastEST {
	
	tag "${chunk_name}"
	publishDir "${params.outdir}/blast_results/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa_ests from fasta_ests 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_ests.collect()
	
	output:
	file blast_result_ests
	
	when:
	params.ESTs != false
		
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa_ests.baseName
	
	"""
	blastn -db $db_name -query $query_fa_ests -max_target_seqs 1 -outfmt 6 -num_threads $params.nthreads > blast_result_ests
	"""
}


/*
 * STEP ESTs.2 - Parse Blast Output
 */

process Blast2QueryTargetEST {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/blast2targets", mode: 'copy'
	
	input:
	file all_blast_results_ests from blast_result_ests.collectFile()
	
	output:
	file query2target_result_uniq_ests into query2target_uniq_result_ests
	
	when:
	params.ESTs != false
	
	script:
	query_tag = ESTs.baseName
	"""
	BlastOutput2QueryTarget.pl $all_blast_results_ests 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq_ests
	"""
}

query2target_uniq_result_ests
	.splitText(by: params.nexonerate, file: true).set{query2target_chunk_ests}	


/*
 * STEP ESTs.3 - Exonerate
 */
 
process RunExonerateEST {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/exonerate/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk_ests
	
	output:
	file 'exonerate.out' into exonerate_result_ests
	
	when:
	params.ESTs != false
	
	script:
	query_tag = ESTs.baseName
	
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_chunk $ESTs $Genome
	"""
}


/*
 * STEP ESTs.4 - Exonerate to Hints
 */
 
process Exonerate2HintsEST {
	
	tag "${query_tag}"
	
	input:
	file exonerate_result_ests
	
	output:
	file exonerate_gff into output_gff_ests
	
	when:
	params.ESTs != false
	
	script:
	query_tag = ESTs.baseName
	
	"""
	grep -v '#' $exonerate_result_ests | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_EST.pl exonerate_gff_lines exonerate_gff
	"""
}

output_gff_ests
 	.collectFile(name: "${params.outdir}/Hints/Hints_ESTs_exonerate.gff")



/*
 * RepeatMasker Block
 */
 
if (params.RM == true) {
	Channel
		.fromPath(Genome)
		.splitFasta(by: params.nrepeats, file: true)
		.set {fasta_rep}
} else {
	fasta_rep = Channel.from(false)
}


/*
 * STEP RepeatMasker.1 - RepeatMasker
 */
 
process RunRepeatMasker {

	tag "${genome_tag}"
	
	publishDir "${params.outdir}/repeatmasker", mode: 'copy'
	
	input:
	file query_fa_rep from fasta_rep 
	
	output:
	file(query_out_rep) into RM_out
	
	script:
	query_out_rep = query_fa_rep + ".out"
	genome_tag = Genome.baseName
	
	"""
	RepeatMasker -species $params.species -pa $params.nthreads $query_fa_rep 
	"""
}


/*
 * STEP RepeatMasker.2 - RepeatMasker - Collect and Clean1
 */
 
process RemoveHeaderRepeatMasker {	
	
	tag "${genome_tag}"
	publishDir "${params.outdir}/repeatmasker", mode: 'copy'
	
	input:
	file "with_header_*" from RM_out.collect()
	
	output:
	file "result_unclean.out" into mergedUNCLEAN

	when:
	params.RM != false
		
	script:
	genome_tag = Genome.baseName
	"""
	tail -n +4 with_header_* > no_header
	cat no_header >> result_unclean.out
	"""
}


/*
 * STEP RepeatMasker.3 - RepeatMasker - Clean2
 */
 
process CleanRepeatMasker {
	
	tag "${genome_tag}"
	
	input:
	file mergedUNCLEAN
	
	output:
	file RepeatMasker_out into RM_2_hints

	when:
	params.RM != false
		
	script:
	genome_tag = Genome.baseName
	
	"""
	grep -v 'with_header' $mergedUNCLEAN | awk 'NF' > RepeatMasker_out
	"""
}


/*
 * STEP RepeatMasker.4 - RepeatMasker to Hints
 */
 
process RepeatMasker2Hints {

	tag "${genome_tag}"
	
	input:
	file RM_2_hints
	
	output:
	file RepeatMasker_hints

	when:
	params.RM != false
	
	script:
	genome_tag = Genome.baseName
		
	"""
	RepeatMasker2hints.pl $RM_2_hints | sort -n -k 1,1 > RepeatMasker_hints
	"""
}

RepeatMasker_hints
	.collectFile(name: "${params.outdir}/Hints/Hints_repeatmasker.gff")



/*
 * RNAseq block
 */


/*
 * Create a channel for input read files
 */
 
 if (params.reads) {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming; read_files_hisat }
} else {
	read_files_fastqc = Channel.from(false)
	read_files_trimming = Channel.from(false)
	read_files_hisat = Channel.from(false)
}


/*
 * STEP RNAseq.1 - FastQC
 */
process runFastqc {
    tag "${prefix}"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
	prefix = reads[0].toString().split("_R1")[0]
    """
    fastqc -q $reads
    """
}


/*
 * STEP RNAseq.2 - Trimgalore
 */
process runTrimgalore {

   tag "${prefix}"
   publishDir "${params.outdir}/trimgalore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (filename.indexOf(".fq") > 0) "$filename"
        }

   input:
   set val(name), file(reads) from read_files_trimming

   output:
   file "*_val_{1,2}.fq" into trimmed_reads
   file "*trimming_report.txt" 
   //into trimgalore_results, trimgalore_logs   
   file "*_fastqc.{zip,html}" 
   //into trimgalore_fastqc_reports
   
   script:
   prefix = reads[0].toString().split("_R1")[0]
   if (params.singleEnd) {
        """
        trim_galore --fastqc --length 36 -q 35 --stringency 1 -e 0.1 $reads
        """
   } else {
        """
        trim_galore --paired --retain_unpaired --fastqc --length 36 -q 35 --stringency 1 -e 0.1 $reads
		"""
   }

}


Channel
	.fromPath(Genome)
	.set { inputMakeHisatdb }

/*
 * STEP RNAseq.3 - Make Hisat2 DB
 */
 
process RunMakeHisatDB {
	
	tag "${prefix}"
	publishDir "${params.outdir}/HisatDB", mode: 'copy'
	
	input:
	file(genome) from inputMakeHisatdb
	
	output:
	file "${dbName}.*.ht2" into hs2_indices
	
	script:
	dbName = genome.baseName
	dbName_1 = dbName + ".1.ht2"
	target = file(dbName_1)
	
	prefix = dbName
    if (!target.exists()) {
		"""
		hisat2-build $genome $dbName -p $params.nthreads
		"""
	}
	
}


/*
 * STEP RNAseq.4 - Hisat2
 */

process RunHisat2 {

	tag "${prefix}"
	publishDir "${params.outdir}/Hisat2", mode: 'copy'
	
	input:
	set val(name), file(reads) from read_files_hisat
	file hs2_indices from hs2_indices.collect()	
	
	output:
	file "*accepted_hits.bam" into accepted_hits2hints, accepted_hits2trinity 
	
	script:
	indexBase = hs2_indices[0].toString() - ~/.\d.ht2/
	ReadsBase = reads[0].toString().split("_R1")[0]
	Read1 = ReadsBase + "_R1.fastq"
	Read2 = ReadsBase + "_R2.fastq"
	prefix = ReadsBase
	
	
	if (params.singleEnd) {
        """
        hisat2 -x $indexBase -U $reads -S alignment_sam -p $params.nthreads
        samtools view -Sb alignment_sam > alignment.bam
        samtools sort alignment.bam > ${prefix}_accepted_hits.bam
        """
   } else {
        """
        hisat2 -x $indexBase -1 $Read1 -2 $Read2 -S alignment_sam -p $params.nthreads
        samtools view -Sb alignment_sam > alignment.bam
        samtools sort alignment.bam > ${prefix}_accepted_hits.bam
		"""
   }
}   


/*
 * STEP RNAseq.5 - Hisat2 into Hints
 */
process Hisat2Hints {

	tag "${prefix}"
	publishDir "${params.outdir}/Hints", mode: 'copy'
	
	input:
	file accepted_hits2hints
	
	output:
	file 'Hints_RNAseq_*.gff'
	
	when:
	params.reads != false
	
	script:
	prefix = accepted_hits2hints[0].toString().split("_accepted")[0]
	
	"""
	bam2hints --intronsonly 0 -p 5 -s 'E' --in=$accepted_hits2hints --out=Hints_RNAseq_${prefix}.gff	
	"""
}


/*
 * STEP RNAseq.6 - Trinity
 */
process runTrinity {

	tag "${prefix}"
	publishDir "${params.outdir}/trinity", mode: 'copy'
	
	input:
	file accepted_hits2trinity
	
	output:
	file "trinity_out_dir/*_trinity.fasta" into trinity_transcripts, trinity_transcripts_2exonerate
	
	when:
	params.reads != false && params.trinity == true
	
	script:
	prefix = accepted_hits2trinity[0].toString().split("_accepted")[0]
	"""
	Trinity --genome_guided_bam $accepted_hits2trinity --genome_guided_max_intron 10000 --CPU $params.nthreads --max_memory 20G
	mv trinity_out_dir/Trinity-GG.fasta trinity_out_dir/${prefix}_trinity.fasta
	"""
}

TrinityChannel = trinity_transcripts.splitFasta(by: params.nblast, file: true)



/*
 * Trinity Transcriptome (Blast + ) Exonerate Block:
 */
 
 
/*
 * STEP RNAseq.7 - Blast
 */
 
process RunBlastTrinity {

	tag "${chunk_name}"
	publishDir "${params.outdir}/blast_trinity/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa from TrinityChannel 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_trinity.collect()
	
	output:
	file blast_result_trinity
	
	when:
	params.reads != false && params.trinity == true	
	
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa.baseName
	
	"""
	blastn -db $db_name -query $query_fa -max_target_seqs 1 -outfmt 6 -num_threads $params.nthreads > blast_result_trinity
	"""
}


/*
 * STEP RNAseq.8 - Parse Blast Output
 */

process BlastTrinity2QueryTarget {
	
	publishDir "${params.outdir}/blast2targets_trinity", mode: 'copy'
	
	input:
	file all_blast_results_trinity from blast_result_trinity.collectFile()
	
	output:
	file query2target_trinity_result_uniq into query2target_trinity_uniq_result
	
	when:
	params.reads != false && params.trinity == true
	
	"""
	BlastOutput2QueryTarget.pl $all_blast_results_trinity 1e-5 query2target_trinity_result
	sort query2target_trinity_result | uniq > query2target_trinity_result_uniq
	"""
} 	

query2target_trinity_uniq_result
	.splitText(by: params.nexonerate, file: true).set{query2target_trinity_chunk}	


/*
 * STEP RNAseq.9 - Exonerate
 */
 
process RunExonerateTrinity {
	
	publishDir "${params.outdir}/exonerate_trinity/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_trinity_chunk from query2target_trinity_chunk
	file trinity_transcripts_2exonerate
	
	output:
	file 'exonerate.out' into exonerate_result_trinity
	
	when:
	params.reads != false && params.trinity == true
		
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_trinity_chunk $trinity_transcripts_2exonerate $Genome
	"""
}


/*
 * STEP RNAseq.10 - Exonerate to Hints
 */
 
process Exonerate2HintsTrinity {
	
	input:
	file exonerate_result_trinity
	
	output:
	file exonerate_trinity_gff into output_trinity_gff, exonerate_trinity_for_hints
	
	when:
	params.reads != false && params.trinity == true	
	
	"""
	grep -v '#' $exonerate_result_trinity | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_trinity.pl exonerate_gff_lines exonerate_trinity_gff
	"""
}

output_trinity_gff
 	.collectFile(name: "${params.outdir}/Hints/Hints_mapped_transcripts.gff")

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    blastn -version > v_blast.txt
    Trinity --version >& v_trinity.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)

/*
 * MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


workflow.onComplete {

    log.info "========================================="
    log.info "Duration:             $workflow.duration"
    log.info "========================================="
        
}
