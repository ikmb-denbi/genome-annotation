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
println x

if ( params.prots ){
	Proteins = file(params.prots)
	x = x + 1
    if( !Proteins.exists() ) exit 1, "Protein file not found: ${Proteins}"
    println "A"
}

if ( params.ESTs ){
	ESTs = file(params.ESTs)
	x = x + 1
    if( !ESTs.exists() ) exit 1, "ESTs file not found: ${ESTs}"
    println "B"
}

if (params.reads){
	x = x + 1
	println "C"
}

println x

if (x == 0) { 
	exit 1, "At least one data file must be especified"
}



println "Hello"




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

if (params.prots) {
	Channel
		.fromPath(Proteins)
		.splitFasta(by: params.nblast, file: true)
		.into {fasta}
} else { 
	fasta = Channel.from(false)
}

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
	
	when:
	params.prots != false
		
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
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
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
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
	
	if (params.qtype == 'protein') {
	"""
	runExonerate_fromBlastHits_prot2genome.pl $hits_chunk $Proteins $Genome
	"""
	} else if (params.qtype == 'EST') {
	"""
	runExonerate_fromBlastHits_est2genome.pl $hits_chunk $Proteins $Genome
	"""
	}
}


/*
 * STEP 5 - Exonerate to Hints
 */
 
process Exonerate2Hints {
	
	tag "${query_tag}"
	
	input:
	file exonerate_result
	
	output:
	file exonerate_gff into output_gff
	
	when:
	params.prots != false
	
	script:
	query_tag = Proteins.baseName
	
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