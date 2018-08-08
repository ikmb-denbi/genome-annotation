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


workflow.onComplete {
        log.info "========================================="
        log.info "Duration:             $workflow.duration"
        log.info "========================================="
}
