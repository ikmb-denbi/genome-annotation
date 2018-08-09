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
	  --trinity						Run transcriptome assembly with trinity and produce hints from the transcripts (true (default) | false)
	  --variant						Specifies whether there are isoforms in the query file ('no_var' (default) | 'var')	
      --nblast						Chunks to divide Blast jobs (default = 10)
      --nexonerate					Chunks to divide Exonerate jobs (default = 10)
	  --nrepeats					Chunks to divide RepeatMasker jobs (default = 2)
	  --species						Species database for RepeatMasker (default = 'mammal')
	  --singleEnd                   Specifies that the input is single end reads (true | false (default))

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
params.variant = "no_var"
params.qtype = "protein"
params.nblast = 10
params.nexonerate = 2
params.nrepeats = 2
params.species = "mammal"
params.name = false
params.singleEnd = false
params.trinity = true


//Script parameters





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



// Create a channel emitting the query fasta file(s), split it in chunks 

if (params.prots) {
	Channel
		.fromPath(Proteins)
		.splitFasta(by: params.nblast, file: true)
		.into {fasta_prots}
} else { 
	fasta_prots = Channel.from(false)
}

//Proteins (Blast + ) Exonerate Block:

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
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	grep -v '#' $exonerate_result_prots | grep 'exonerate:protein2genome:local' > exonerate_gff_lines
	Exonerate2GFF_protein.pl exonerate_gff_lines $params.variant exonerate_gff
	"""
}

output_gff_prots
 	.collectFile(name: "${params.outdir}/Hints_proteins_exonerate.gff")


if (params.ESTs) {
	Channel
		.fromPath(ESTs)
		.splitFasta(by: params.nblast, file: true)
		.into {fasta_ests}
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
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
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
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_chunk_ests $Proteins $Genome
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
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	"""
	grep -v '#' $exonerate_result_ests | grep 'exonerate:est2genome' > exonerate_gff_lines
	Exonerate2GFF_EST.pl exonerate_gff_lines $params.variant exonerate_gff
	"""
}

output_gff_ests
 	.collectFile(name: "${params.outdir}/Hints_ESTs_exonerate.gff")