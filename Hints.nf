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
      --query						Proteins from other species
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. docker / aws

    Options:
	  --variant						Specifies whether there are isoforms in the query file ('no_var' (default) | 'var')
      --qtype						Query type: ('protein' (default) | 'EST')
      --nblast						Chunks to divide Blast jobs (default = 10)
      --nexonerate					Chunks to divide Exonerate jobs (default = 10)
	  --nrepeats					Chunks to divide RepeatMasker jobs (default = 2)
	  --species						Species database for RepeatMasker (default = 'mammal')
	  --singleEnd                   Specifies that the input is single end reads

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
params.variant = "no_var"
params.qtype = "protein"
params.nblast = 10
params.nexonerate = 10
params.nrepeats = 2
params.species = "mammal"
params.name = false
params.outdir = "./Hints_output"



// Validate inputs
if ( params.genome ){
    fasta = file(params.genome)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.genome}"
}

if ( params.query ){
    fasta = file(params.query)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.query}"
}


//Script parameters
Queries = file(params.query)
Genome = file(params.genome)


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
summary['Query']		= Queries
summary['Query type']	= params.qtype
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
	set file(db_nhr),file(db_nin),file(db_nsq) into blast_db, blast_db_trinity
	
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



// Create a channel emitting the query fasta file(s), split it in chunks 

Channel
	.fromPath(Queries)
	.splitFasta(by: params.nblast, file: true)
	.into {fasta}


//Proteins (Blast + ) Exonerate Block:

/*
 * STEP 2 - Blast
 */
 
process RunBlast {
	
	tag "${chunk_name}"
	publishDir "${params.outdir}/blast_results/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa from fasta 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db.collect()
	
	output:
	file blast_result
		
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa.baseName
	
	if (params.qtype == 'protein') {
		"""
		tblastn -db $db_name -query $query_fa -max_target_seqs 1 -outfmt 6 > blast_result
		"""
	} else if (params.qtype == 'EST') {
		"""
		blastn -db $db_name -query $query_fa -max_target_seqs 1 -outfmt 6 > blast_result
		"""
	}
}



/*
 * STEP 3 - Parse Blast Output
 */

process Blast2QueryTarget {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/blast2targets", mode: 'copy'
	
	input:
	file all_blast_results from blast_result.collectFile()
	
	output:
	file query2target_result_uniq into query2target_uniq_result
	
	script:
	query_tag = Queries.baseName
	"""
	BlastOutput2QueryTarget.pl $all_blast_results 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq
	"""
}

query2target_uniq_result
	.splitText(by: params.nexonerate, file: true).set{query2target_chunk}	


/*
 * STEP 4 - Exonerate
 */
 
process RunExonerate {
	
	tag "${query_tag}"
	publishDir "${params.outdir}/exonerate/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk
	
	output:
	file 'exonerate.out' into exonerate_result
	
	script:
	query_num = hits_chunk.toString().split(".")[-1]
	query_tag = Queries.baseName + " " + query_num
	
	
	if (params.qtype == 'protein') {
	"""
	runExonerate_fromBlastHits_prot2genome.pl $hits_chunk $Queries $Genome
	"""
	} else if (params.qtype == 'EST') {
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_chunk $Queries $Genome
	"""
	}
}


/*
 * STEP 5 - Exonerate to Hints
 */
 
process Exonerate2Hints {
	
	input:
	file exonerate_result
	
	output:
	file exonerate_gff into output_gff
	
	script:
	if (params.qtype == 'protein') {
	"""
	grep -v '#' $exonerate_result | grep 'exonerate:protein2genome:local' > exonerate_gff_lines
	Exonerate2GFF_protein.pl exonerate_gff_lines $params.variant exonerate_gff
	"""
	} else if (params.qtype == 'EST') {
	"""
	grep -v '#' $exonerate_result | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_EST.pl exonerate_gff_lines $params.variant exonerate_gff
	"""
	}
}

output_gff
 	.collectFile(name: "${params.outdir}/Hints_${params.qtype}_exonerate.gff")


/*
 * STEP 6 - GenomeThreader
 */
 
process RunGenomeThreader {
	
	publishDir "${params.outdir}/genomethreader", mode: 'copy'
		
	output:
	file output_gth
	
	script:
	if (params.qtype == 'protein') {
	"""
	gth -genomic $Genome -protein $Queries -gff3out -intermediate -o output_gth
	"""
	} else if (params.qtype == 'EST') {
	"""
	gth -genomic $Genome -cdna $Queries -gff3out -intermediate -o output_gth
	"""
	}
}


/*
 * STEP 7 - GenomeThreader to Hints
 */
 
process GenomeThreader2Hints {
	
	input:
	file not_clean_gth from output_gth
	
	output:
	file gth_hints
	
	"""
	gt gff3 -addintrons yes -setsource gth -tidy yes -addids no $not_clean_gth > not_clean_gth_wIntrons
	grep -v '#' not_clean_gth_wIntrons > no_hash_gth
	GTH_rename_splitoutput.pl no_hash_gth > clean_gth
	grep -e 'CDS' -e 'exon' -e 'intron' clean_gth | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > gth_hints
	"""
}

gth_hints
	.collectFile(name: "${params.outdir}/Hints_${params.qtype}_genomethreader.gff")


//RepeatMasker Block

Channel
	.fromPath(Genome)
	.splitFasta(by: params.nrepeats, file: true)
	.set {fasta_rep}


/*
 * STEP 8 - RepeatMasker
 */
 
process RunRepeatMasker {

	publishDir "${params.outdir}/repeatmasker", mode: 'copy'
	
	input:
	file query_fa_rep from fasta_rep 
	
	output:
	file(query_out_rep) into RM_out
	
	script:
	query_out_rep = query_fa_rep + ".out"
	
	"""
	RepeatMasker -species $params.species $query_fa_rep
	"""
}

/*
 * STEP 9 - RepeatMasker - Collect and Clean1
 */
 
process RemoveHeaderRepeatMasker {	
	
	publishDir "${params.outdir}/repeatmasker", mode: 'copy'
	
	input:
	file "with_header_*" from RM_out.collect()
	
	output:
	file "result_unclean.out" into mergedUNCLEAN
	
	"""
	tail -n +4 with_header_* > no_header
	cat no_header >> result_unclean.out
	"""
}


/*
 * STEP 10 - RepeatMasker - Clean2
 */
 
process CleanRepeatMasker {
	
	input:
	file mergedUNCLEAN
	
	output:
	file RepeatMasker_out into RM_2_hints
	"""
	grep -v 'with_header' $mergedUNCLEAN | awk 'NF' > RepeatMasker_out
	"""
}


/*
 * STEP 11 - RepeatMasker to Hints
 */
 
process RepeatMasker2Hints {

	input:
	file RM_2_hints
	
	output:
	file RepeatMasker_hints
	
	"""
	RepeatMasker2hints.pl $RM_2_hints | sort -n -k 1,1 > RepeatMasker_hints
	"""
}

RepeatMasker_hints
	.collectFile(name: "${params.outdir}/Hints_repeatmasker.gff")



/*
 * RNAseq block
 */
 
 
/*
 * Create a channel for input read files
 */
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming; read_files_hisat }

/*
 * STEP 12 - FastQC
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
 * STEP 13 - Trimgalore
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
 * STEP 14 - Make Hisat2 DB
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
		hisat2-build $genome $dbName
		"""
	}
	
}


/*
 * STEP 15 - Hisat2
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
        hisat2 -x $indexBase -U $reads -S alignment_sam
        samtools view -Sb alignment_sam > alignment.bam
        samtools sort alignment.bam > ${prefix}_accepted_hits.bam
        """
   } else {
        """
        hisat2 -x $indexBase -1 $Read1 -2 $Read2 -S alignment_sam
        samtools view -Sb alignment_sam > alignment.bam
        samtools sort alignment.bam > ${prefix}_accepted_hits.bam
		"""
   }
}   

/*
 * STEP 16 - Hisat2 into Hints
 */
process Hisat2Hints {

	tag "${prefix}"
	publishDir "${params.outdir}", mode: 'copy'
	
	input:
	file accepted_hits2hints
	
	output:
	file 'Hints_RNAseq_*.gff'
	
	script:
	prefix = accepted_hits2hints[0].toString().split("_accepted")[0]
	
	"""
	bam2hints --intronsonly 0 -p 5 -s 'E' --in=$accepted_hits2hints --out=Hints_RNAseq_${prefix}.gff	
	"""
}

/*
 * STEP 17 - Trinity
 */
process runTrinity {

	tag "${prefix}"
	publishDir "${params.outdir}/trinity", mode: 'copy'
	
	input:
	file accepted_hits2trinity
	
	output:
	file "trinity_out_dir/*_trinity.fasta" into trinity_transcripts, trinity_transcripts_2exonerate
	
	script:
	prefix = accepted_hits2trinity[0].toString().split("_accepted")[0]
	"""
	Trinity --genome_guided_bam $accepted_hits2trinity --genome_guided_max_intron 10000 --CPU 1 --max_memory 5G
	mv trinity_out_dir/Trinity-GG.fasta trinity_out_dir/${prefix}_trinity.fasta
	"""
}

TrinityChannel = trinity_transcripts.splitFasta(by: params.nblast, file: true)


//Proteins (Blast + ) Exonerate Block:

/*
 * STEP 18 - Blast
 */
 
process RunBlastTrinity {

	publishDir "${params.outdir}/blast_trinity/${chunk_name}", mode: 'copy'
	
	input:
	file query_fa from TrinityChannel 
	set file(blastdb_nhr),file(blast_nin),file(blast_nsq) from blast_db_trinity.collect()
	
	output:
	file blast_result_trinity
		
	script: 

	db_name = blastdb_nhr.baseName
	chunk_name = query_fa.baseName
	
	"""
	blastn -db $db_name -query $query_fa -max_target_seqs 1 -outfmt 6 > blast_result_trinity
	"""
}

/*
 * STEP 19 - Parse Blast Output
 */

process BlastTrinity2QueryTarget {
	
	publishDir "${params.outdir}/blast2targets_trinity", mode: 'copy'
	
	input:
	file all_blast_results_trinity from blast_result_trinity.collectFile()
	
	output:
	file query2target_trinity_result_uniq into query2target_trinity_uniq_result
	
	"""
	BlastOutput2QueryTarget.pl $all_blast_results_trinity 1e-5 query2target_trinity_result
	sort query2target_trinity_result | uniq > query2target_trinity_result_uniq
	"""
} 	

query2target_trinity_uniq_result
	.splitText(by: params.nexonerate, file: true).set{query2target_trinity_chunk}	


/*
 * STEP 20 - Exonerate
 */
 
process RunExonerateTrinity {

	publishDir "${params.outdir}/exonerate_trinity/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_trinity_chunk from query2target_trinity_chunk
	file trinity_transcripts_2exonerate
	
	output:
	file 'exonerate.out' into exonerate_result_trinity
	
	script:
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_trinity_chunk $trinity_transcripts_2exonerate $Genome
	"""
}


/*
 * STEP 21 - Exonerate to Hints
 */
 
process Exonerate2HintsTrinity {
	
	input:
	file exonerate_result_trinity
	
	output:
	file exonerate_trinity_gff into output_trinity_gff, exonerate_trinity_for_hints
	
	script:
	"""
	grep -v '#' $exonerate_result_trinity | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_trinity.pl exonerate_gff_lines exonerate_trinity_gff
	"""
}

output_trinity_gff
 	.collectFile(name: "${params.outdir}/Hints_mapped_transcripts.gff")


workflow.onComplete {
        log.info "========================================="
        log.info "Duration:             $workflow.duration"
        log.info "========================================="
}
