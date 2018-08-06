#!/usr/bin/env nextflow
/*
========================================================================================
                         NF-hints
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

    nextflow run NF-hints --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Genome reference
      --queries						Proteins from other species
      --qtype						Query type: protein/EST
      --nblast						Chunks to divide Blast jobs
      --nexonerate					Chunks to divide Exonerate jobs
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
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

// Configurable variables
params.name = false
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
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
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
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
if(params.email) summary['E-mail Address'] = params.email
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
	.fromPath(params.genome)
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
	set file(db_nhr),file(db_nin),file(db_nsq) into blast_db
	
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
	.fromPath(params.queries)
	.splitFasta(by: params.nblast, file: true)
	.ifEmpty { exit 1, "Could not find proteins file" }
	.into {fasta}


//Proteins (Blast + ) Exonerate Block:

/*
 * STEP 2 - Blast
 */
 
process RunBlast {

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
	
	publishDir "${params.outdir}/blast2targets", mode: 'copy'
	
	input:
	file all_blast_results from blast_result.collectFile()
	
	output:
	file query2target_result_uniq into query2target_uniq_out, query2target_uniq_result
	
	"""
	BlastOutput2QueryTarget.pl $all_blast_results 1e-5 query2target_result
	sort query2target_result | uniq > query2target_result_uniq
	"""
}

query2target_uniq_out
	.collectFile(name: "${params.outdir}/Blast_output.txt") 	

query2target_uniq_result
	.splitText(by: params.nexonerate, file: true).set{query2target_chunk}	


/*
 * STEP 5 - Exonerate
 */
 
process RunExonerate {

	publishDir "${params.outdir}/exonerate/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk
	file Genome from file(params.genome)
	file Queries from file(params.queries)
	
	output:
	file 'exonerate.out' into exonerate_result
	
	script:
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
 * STEP 6 - Exonerate to GFF
 */
 
process Exonerate2Gff {
	
	input:
	file exonerate_result
	
	output:
	file exonerate_gff into output_gff, exonerate_for_hints
	
	script:
	if (params.qtype == 'protein') {
	"""
	grep -v '#' $exonerate_result | grep 'exonerate:protein2genome:local' > exonerate_gff_lines
	perl /ifs/data/nfs_share/sukmb415/scripts/Exonerate2GFF_protein.pl exonerate_gff_lines $VARIANT exonerate_gff
	"""
	} else if (params.qtype == 'EST') {
	"""
	grep -v '#' $exonerate_result | grep 'exonerate:est2genome' > exonerate_gff_lines
	perl /ifs/data/nfs_share/sukmb415/scripts/Exonerate2GFF_EST.pl exonerate_gff_lines $VARIANT exonerate_gff
	"""
	}
}

output_gff
 	.collectFile(name: "${params.outdir}/Exonerate_output.gff")





/*
 * Create a channel for input read files
 */
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming }



/*
 * STEP X - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

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
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP X - MultiQC
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



/*
 * STEP X - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



workflow.onComplete {
        log.info "========================================="
        log.info "Duration:             $workflow.duration"
        log.info "========================================="
}
