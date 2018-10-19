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
----------------------------------------------------------------------------------------
*/


def helpMessage() {
	log.info"""
  =================================================================
   IKMB - de.NBI | Genome Annotation | Pipeline v${params.version}
  =================================================================
  Usage:

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
    --ninterpro		Chunks (# of sequences) to divide InterPro jobs [ default = 50 ]
    --nthreads		Number of cpus for programs that allow multi-threaded mode [default = 1]	

    Other options:
    --singleEnd		Specifies that the input is single end reads [ true | false (default) ]
    --outdir		The output directory where the results will be saved [ default = 'Hints_augustus_output' ]
    --allHints		Name of final GFF file with all hints [ default = 'AllHints.gff' ]
    --addHints		Additional hints file (in GFF format), to be concatenated to the resulting hints before running augustus [ default = 'false' ]
    -name			Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
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



//Default variables:
params.prots = false
params.ESTs = false
params.reads = false

params.trinity = true
params.gth = true
params.RM = true
params.augustus = true
params.funAnnot = true

params.species = "mammal"
params.model = "human"
params.UTR = 'off'
params.isof = 'false'
params.augCfg = false
params.uniprot = false

params.nblast = 20
params.nexonerate = 5
params.nrepeats = 1
params.ninterpro = 10
params.nthreads = 1

params.singleEnd = false
params.outdir = "Hints_augustus_output"
params.allHints = "AllHints.gff"
params.addHints = false
params.name = false


AllHints = file(params.allHints)

GFF3_RUBYscript = file(workflow.projectDir + "/bin/augustus_add_exons.rb")
ADDANNO_RIBYscript = file(workflow.projectDir + "/bin/gff_add_annie_functions.rb")
CUR_DIR = "$PWD"

if(params.augCfg == false) {
	AUG_CONF = file(workflow.projectDir + "/bin/augustus_default.cfg")
} else {
	AUG_CONF = params.augCfg
}

if(params.uniprot == false) {
	UNIPROTDB = file(workflow.projectDir + "/bin/Eumetazoa_UniProt_reviewed_evidence.fa")
} else {
	UNIPRTOTDB = params.uniprot
}


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

// Is there already a Hints file from a previous run?
if(AllHints.exists() ) {
	exit 1, "$AllHints already exists, please remove it or give a different name with --allHints"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
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
summary['Pipeline Name']  = 'NF-Hints-Augustus'
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
	
// Check if the blast db already exists - if not, we create it

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
	trigger_prot_exonerate = Channel.create()
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
	tblastn -db $db_name -query $query_fa_prots -max_target_seqs 1 -outfmt 6 > blast_result_prots
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
	file 'prot_exonerate_hints.done' into trigger_prot_exonerate
	
	script:
	query_tag = Proteins.baseName
	
	if (params.prots != false) {
	"""
	grep -v '#' $exonerate_result_prots | grep 'exonerate:protein2genome:local' > exonerate_gff_lines
	Exonerate2GFF_protein.pl exonerate_gff_lines exonerate_gff
	cat exonerate_gff >> $AllHints
	touch prot_exonerate_hints.done
	"""
	}
}

output_gff_prots
 	.collectFile(name: "${params.outdir}/Hints/Hints_proteins_exonerate.gff")


if (params.gth == false) {
	trigger_prot_gth = Channel.create()
}

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
	file 'prot_gth_hints.done' into trigger_prot_gth
	
	script:
	query_tag = Proteins.baseName
	
	if (params.prots != false && params.gth != false) {
	"""
	gt gff3 -addintrons yes -setsource gth -tidy yes -addids no $not_clean_gth > not_clean_gth_wIntrons
	grep -v '#' not_clean_gth_wIntrons > no_hash_gth
	GTH_rename_splitoutput.pl no_hash_gth > clean_gth
	grep -e 'CDS' -e 'exon' -e 'intron' clean_gth | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > gth_hints
	cat gth_hints >> $AllHints
	touch prot_gth_hints.done
	"""
	} 
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
	trigger_est_exonerate = Channel.create()
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
	blastn -db $db_name -query $query_fa_ests -max_target_seqs 1 -outfmt 6 > blast_result_ests
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
	file 'est_exonerate_hints.done' into trigger_est_exonerate
	
	script:
	query_tag = ESTs.baseName
	
	if (params.ESTs != false) {
	"""
	grep -v '#' $exonerate_result_ests | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_EST.pl exonerate_gff_lines exonerate_gff
	cat exonerate_gff >> $AllHints
	touch est_exonerate_hints.done
	"""
	}
}

output_gff_ests
 	.collectFile(name: "${params.outdir}/Hints/Hints_ESTs_exonerate.gff")



/*
 * RepeatMasker Block
 */
 
if (params.RM != false) {
	Channel
		.fromPath(Genome)
		.splitFasta(by: params.nrepeats, file: true)
		.set {fasta_rep}
} else {
	fasta_rep = Channel.from(false)
	trigger_RM = Channel.create()
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
	
	when:
	params.RM != false

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
	file RepeatMasker_gff into RepeatMasker_hints, RepeatMasker_2concatenate
	file 'RM_hints.done' into trigger_RM
	
	script:
	genome_tag = Genome.baseName
	
	if (params.RM != false) {	
	"""
	RepeatMasker2hints.pl $RM_2_hints | sort -n -k 1,1 > RepeatMasker_gff
	cat RepeatMasker_gff >> $AllHints
	touch RM_hints.done
	"""
	}
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
         .into { read_files_fastqc; read_files_trimming }
} else {
	read_files_fastqc = Channel.from(false)
	read_files_trimming = Channel.from(false)
	trigger_RNAseq = Channel.create()
	
}


/*
 * STEP RNAseq.1 - FastQC
 */
process RunFastqc {
    tag "${prefix}"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" 
    //into fastqc_results

	when:
	params.reads != false

    script:
	prefix = reads[0].toString().split("_R1")[0]
    """
    fastqc -q $reads
    """
}


/*
 * STEP RNAseq.2 - Trimgalore
 */
process RunTrimgalore {

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
  // file "*trimming_report.txt" into trimgalore_results, trimgalore_logs   
   file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
 	
	when:
	params.reads != false
  
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
	
	when:
	params.reads != false

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
	file reads from trimmed_reads
	file hs2_indices from hs2_indices.collect()	
	
	output:
	file "*accepted_hits.bam" into accepted_hits2hints, accepted_hits2trinity 
	
	when:
	params.reads != false
	
	script:
	indexBase = hs2_indices[0].toString() - ~/.\d.ht2/
	ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

	prefix = ReadsBase
	
	
	if (params.singleEnd) {
        """
        hisat2 -x $indexBase -U $reads -S alignment_sam -p $params.nthreads
        samtools view -Sb alignment_sam > alignment.bam
        samtools sort alignment.bam > ${prefix}_accepted_hits.bam
        """
   } else {
        """
        hisat2 -x $indexBase -1 ${reads[0]} -2 ${reads[1]} -S alignment_sam -p $params.nthreads
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
	file 'RNAseq_hints.done' into trigger_RNAseq
	
	script:
	prefix = accepted_hits2hints[0].toString().split("_accepted")[0]
	
	if (params.reads != false) {
	"""
	bam2hints --intronsonly 0 -p 5 -s 'E' --in=$accepted_hits2hints --out=Hints_RNAseq_${prefix}.gff	
	cat Hints_RNAseq_${prefix}.gff >> $AllHints
	touch RNAseq_hints.done
	"""
	}
}

if (params.trinity == false) {
	trigger_trinity = Channel.create()
}

/*
 * STEP RNAseq.6 - Trinity
 */
process RunTrinity {

	publishDir "${params.outdir}/trinity", mode: 'copy'
	
	input:
	file hisathits from accepted_hits2trinity.collect()
	
	output:
	file "trinity_out_dir/Transcripts_trinity.fasta" into trinity_transcripts, trinity_transcripts_2exonerate
	
	when:
	params.reads != false && params.trinity == true
	
	script:
	
	"""
	samtools merge merged.bam $hisathits
	samtools sort merged.bam > sorted.bam
	Trinity --genome_guided_bam sorted.bam --genome_guided_max_intron 10000 --CPU $params.nthreads --max_memory 20G
	mv trinity_out_dir/Trinity-GG.fasta trinity_out_dir/Transcripts_trinity.fasta
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
	
	publishDir "${params.outdir}/blast_trinity/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa from TrinityChannel 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_trinity.collect()
	
	output:
	file blast_trinity
	
	when:
	params.reads != false && params.trinity == true	
	
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa.baseName
	
	"""
	blastn -db $db_name -query $query_fa -max_target_seqs 1 -outfmt 6 > blast_trinity
	"""
}


/*
 * STEP RNAseq.8 - Parse Blast Output
 */

process BlastTrinity2QueryTarget {
	
	publishDir "${params.outdir}/blast2targets_trinity", mode: 'copy'
	
	input:
	file all_blast_results_trinity from blast_trinity.collectFile()
	
	output:
	file query2target_trinity_uniq into query2target_trinity_uniq_result
	
	when:
	params.reads != false && params.trinity == true
	
	script:

	"""
	BlastOutput2QueryTarget.pl $all_blast_results_trinity 1e-5 query2target_trinity_result
	sort query2target_trinity_result | uniq > query2target_trinity_uniq
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
	
	script:
		
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
	file Hints_trinity_gff into Hints_trinity2concatenate, Hints_trinity_mapped_gff
	file 'trinity_hints.done' into trigger_trinity
	
	script:
	if (params.reads != false && params.trinity == true) {
	"""
	grep -v '#' $exonerate_result_trinity | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_trinity.pl exonerate_gff_lines Hints_trinity_gff
	cat Hints_trinity_gff >> $AllHints
	touch trinity_hints.done
	"""
	}
}

Hints_trinity_mapped_gff
	.collectFile(name: "${params.outdir}/Hints/Hints_trinity_mapped.gff")


/*
 * Augustus Block:
 */
 
 		
/*
 * STEP Augustus.1 - Genome Annotation
 */
process RunAugustus {

	input:
	file a from trigger_prot_exonerate.ifEmpty()
	file b from trigger_prot_gth.ifEmpty()
	file c from trigger_est_exonerate.ifEmpty()
	file d from trigger_RM.ifEmpty()
	file e from trigger_RNAseq.ifEmpty()
	file f from trigger_trinity.ifEmpty()	
	
	output:
	file Augustus_out into augustus_out_gff, augustus_2gff3, augustus_2prots
	
	when:
	params.augustus != false
	
	script:
	if (params.addHints == false) {
	"""
	augustus --species=$params.model --UTR=$params.UTR --alternatives-from-evidence=$params.isof --extrinsicCfgFile=$AUG_CONF --hintsfile=$AllHints $Genome > Augustus_out
	"""
	} else {
	AdditionalHints = "$CUR_DIR" + "/" +  "$params.addHints"
	"""
	cat $AllHints $AdditionalHints >> combinedHints
	augustus --species=$params.model --UTR=$params.UTR --alternatives-from-evidence=$params.isof --extrinsicCfgFile=$AUG_CONF --hintsfile=combinedHints $Genome > Augustus_out
	"""
	}
}

augustus_out_gff
	.collectFile( name: "${params.outdir}/Augustus_out.gff" )


/*
 * STEP Augustus.2 - Get GFF3 file
 */

process Augustus2Gff3 {
	
	input:
	file augustus2parse from augustus_2gff3
	
	output:
	file augustus_gff3 into augustus_gff3_out, augustus_gff32annie, augustus_gff32interpro, augustus_gff3annotate
	
	when:
	params.augustus != false
	
	"""
	grep -v '#' $augustus2parse | sed 's/transcript/mRNA/' > augustus_clean
	ruby $GFF3_RUBYscript -i augustus_clean > augustus_gff3
	"""
}

augustus_gff3_out
	.collectFile( name: "${params.outdir}/Augustus.gff3" )

/*
 * STEP Augustus.3 - Get Protein Sequences
 */

process Augustus2Proteins {
	
	input:
	file augustus2parse from augustus_2prots
	
	output:
	file '*.aa' into augustus_proteins, augustus_prots2annie, augustus_prots2interpro
	
	when:
	params.augustus != false
	
	"""
	getAnnoFasta.pl $augustus2parse
	"""
}

augustus_proteins
	.collectFile( name: "${params.outdir}/Augustus_proteins.fa" )


/*
 * Functional Annotation Block:
 */


Channel
	.fromPath(UNIPROTDB)
	.set { inputMakeblastdb }
	 
/*
 * STEP Functional Annotation.1 - Uniprot BlastDB
 */
 
// We check if the blast db already exists - if not, we create it

process RunMakeBlastDBFunAnno {
	
	publishDir "${params.outdir}/FunAnnoBlastDB", mode: 'copy'
	
	input:
	file(uniprot_fa) from inputMakeblastdb
	
	output:
	set file(db_phr),file(db_pin),file(db_psq) into blast_db
	
	when:
	params.funAnnot != false
	
	script:
	dbName = uniprot_fa.baseName
	db_phr = dbName + ".phr"
	db_pin = dbName + ".pin"
	db_psq = dbName + ".psq"

	target = file(db_phr)
	
    if (!target.exists()) {
		"""
		makeblastdb -in $uniprot_fa -dbtype prot -out $dbName
		"""
	}
	
}


/*
 * STEP Functional Annotation.2 - BlastP of annotated proteins against UniProtDB
 */

augustus_prots2annie
	.splitFasta(by: params.nblast, file: true).set{proteinChunkBlast}

process RunBlastpFunAnno {

	publishDir "${params.outdir}/blast_results_FunnAnno/${chunk_name}", mode: 'copy'
	input:
    file fasta from proteinChunkBlast
	set file(blastdb_phr),file(blastdb_pin),file(blastdb_psq) from blast_db.collect()

	output:
    file blast_result_funAnno

	when:
	full == true
	params.funAnnot != false

	script:
	db_name = blastdb_phr.baseName
	chunk_name = fasta.baseName
 
	"""
	blastp -query $fasta -db $db_name -evalue 0.01 -outfmt 6 -num_threads 4 > blast_result_funAnno
	"""

}

mergedBlast = blast_result_funAnno.collectFile(name:'mergedBlast')

/*
 * STEP Functional Annotation.3 - Annie on BlastP results
 */
 
process RunAnnieBlast {

	publishDir "${params.outdir}/annie"

	input:
	file blast from mergedBlast
	file augustus_gff32annie

	output:
	file annie_report into outputAnnieBlast
	
	when:
	params.funAnnot != false
	
	script:
	annie_report = "blastp.annie"

	"""
	annie.py -db $UNIPROTDB -b $blast -g $augustus_gff32annie -o $annie_report
	"""

}

/*
 * STEP Functional Annotation.4 - InterPro
 */
 
augustus_prots2interpro
	.splitFasta(by: params.ninterpro, file: true).set{proteinChunkInterpro}

process RunInterproscan {
   
	publishDir "${params.outdir}/interpro_results/${chunk_name}", mode: 'copy'
   
	input:
	file fasta from proteinChunkInterpro

	output:
	file interpro into outputInterpro
  
	when:
	full == true
	params.funAnnot != false
	
	script:
	interpro = "interpro.tsv"
	chunk_name = fasta.baseName

	"""
	interproscan.sh -appl Pfam -i $fasta -b interpro -iprlookup -goterms -pa -dp -f tsv 
	"""

}

mergedInterPro = outputInterpro.collectFile(name:'mergedInterPro')

/*
 * STEP Functional Annotation.5 - Annie on InterPro results
 */
 
process RunAnnieInterpro {

	publishDir "${params.outdir}/annie"

	input:
	file ipr from mergedInterPro
	file augustus_gff32interpro
      
	output:
	file annie_report into outputAnnieInterpro

	when:
	params.funAnnot != false
	
	script:
	annie_report = "interpro.annie"

	"""
	annie.py -ipr $ipr -o $annie_report -g $augustus_gff32interpro
	"""

}

/*
 * STEP Functional Annotation.6 - Functional annotation to GFF3
 */
 
process RunFunctionsToGFF {
 
	input:
	file interpro from outputAnnieInterpro
	file blast from outputAnnieBlast
	file augustus_gff3annotate
     
	output:
	file annotated_gff 
	
	when:
	params.funAnnot != false
	
  	"""
	cat $interpro $blast > annie.txt
	ruby $ADDANNO_RIBYscript -g $augustus_gff3annotate -a annie.txt > annotated_gff
  	"""

}

annotated_gff
	.collectFile( name: "${params.outdir}/Augustus_withFunctions.gff3" )

/*
 * Parse software version numbers
 */
 
process Get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
  
    if (params.gth == true)
    	"""
    	fastqc --version > v_fastqc.txt
    	trim_galore --version &> v_trim_galore.txt
    	hisat2 --version > v_hisat2.txt
    	blastn -version > v_blast.txt
    	#exonerate -v > v_exonerate.txt
    	gth -version > v_gth.txt
    	RepeatMasker -v > v_rm.txt
    	#Trinity --version > v_trinity.txt
    	echo $params.version > v_pipeline.txt
    	echo $workflow.nextflow.version > v_nextflow.txt    
    	multiqc --version > v_multiqc.txt
    	scrape_software_versions.py > software_versions_mqc.yaml
    	"""
    else
    	"""
    	fastqc --version > v_fastqc.txt
    	trim_galore --version &> v_trim_galore.txt
    	hisat2 --version > v_hisat2.txt
    	blastn -version > v_blast.txt
    	#exonerate -v > v_exonerate.txt
    	RepeatMasker -v > v_rm.txt
    	#Trinity --version > v_trinity.txt
    	echo $params.version > v_pipeline.txt
    	echo $workflow.nextflow.version > v_nextflow.txt    
    	multiqc --version > v_multiqc.txt
    	scrape_software_versions.py > software_versions_mqc.yaml
   	 	"""
}


params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
multiqc_config = file(params.multiqc_config)

/*
 * MultiQC
 */
process Multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
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
