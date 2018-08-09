#!/usr/bin/env nextflow
/*
========================================================================================
                         NF-RNAseq
========================================================================================
 NF-RNAseq Analysis Pipeline. Started 2018-08-03.
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

    nextflow run NF-hints --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Genome reference
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    Other options:
      --outdir                      The output directory where the results will be saved
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

//Script parameters
Queries = file(params.query)
Genome = file(params.genome)
params.nblast = 1
params.nexonerate = 2

Channel
	.fromPath(Genome)
	.set { inputMakeblastdb }
	
// We check if the blast db already exists - if not, we create it

/*
 * STEP 1 - Make Blast DB
 */
 
process RunMakeBlastDB {
	
	publishDir "${params.outdir}/BlastDB", mode: 'copy'
	
	input:
	file(genome) from inputMakeblastdb
	
	output:
	set file(db_nhr),file(db_nin),file(db_nsq) into blast_db_trinity
	
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
 * Create a channel for input read files
 */
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming; read_files_hisat }

/*
 * STEP 1 - FastQC
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
 * STEP 2 - Trimgalore
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
 * STEP 3 - Make Hisat2 DB
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
 * STEP 4 - Hisat2
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
	prefix = ReadsBase + "_vs_" + indexBase
	
	
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
 * STEP 5 - Bam into Hints
 */
process Bam2Hints {

	tag "${prefix}"
	publishDir "${params.outdir}", mode: 'copy'
	
	input:
	file accepted_hits2hints
	
	output:
	file '*_hints.gff'
	
	script:
	prefix = accepted_hits2hints[0].toString().split("_accepted")[0]
	
	"""
	bam2hints --intronsonly 0 -p 5 -s 'E' --in=$accepted_hits2hints --out=${prefix}_hints.gff	
	"""
}

/*
 * STEP 6 - Trinity
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
 * STEP 7 - Blast
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
 * STEP 8 - Parse Blast Output
 */

process BlastTrinity2QueryTarget {
	
	publishDir "${params.outdir}/blast2targets_trinity", mode: 'copy'
	
	input:
	file all_blast_results_trinity from blast_result_trinity.collectFile()
	
	output:
	file query2target_trinity_result_uniq into query2target_trinity_uniq_out, query2target_trinity_uniq_result
	
	"""
	BlastOutput2QueryTarget.pl $all_blast_results_trinity 1e-5 query2target_trinity_result
	sort query2target_trinity_result | uniq > query2target_trinity_result_uniq
	"""
}

query2target_trinity_uniq_out
	.collectFile(name: "${params.outdir}/Blast_trinity_output.txt") 	

query2target_trinity_uniq_result
	.splitText(by: params.nexonerate, file: true).set{query2target_trinity_chunk}	


/*
 * STEP 9 - Exonerate
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
 * STEP 10 - Exonerate to GFF
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
 	.collectFile(name: "${params.outdir}/Exonerate_trinity_transcript_hints.gff")


workflow.onComplete {
        log.info "========================================="
        log.info "Duration:             $workflow.duration"
        log.info "========================================="
}
