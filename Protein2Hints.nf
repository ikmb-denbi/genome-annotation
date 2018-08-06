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

    nextflow run Protein2Hints.nf --genome 'Genome.fasta' --query 'Proteins.fasta' -profile docker

    Mandatory arguments:
      --genome                      Genome reference
      --query						Proteins from other species
      -profile                      Hardware config to use. docker / aws

    Options:
	  --variant						Specifies whether there are isoforms in the query file ('no_var' (default) | 'var')
      --qtype						Query type: ('protein' (default) | 'EST')
      --nblast						Chunks to divide Blast jobs (default = 10)
      --nexonerate					Chunks to divide Exonerate jobs (default = 10)
	  --nrepeats					Chunks to divide RepeatMasker jobs (default = 2)
	  --species						Species database for RepeatMasker (default = 'mammal'(

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
params.outdie = "Protein2Hints_output"

//Script parameters
Queries = file(params.query)
Genome = file(params.genome)


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
summary['Fasta Ref']    = params.fasta
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
	.fromPath(Queries)
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
 * STEP 4 - Exonerate
 */
 
process RunExonerate {

	publishDir "${params.outdir}/exonerate/${hits_chunk}", mode: 'copy'
	
	input:
	file hits_chunk from query2target_chunk
	
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
 * STEP 5 - Exonerate to GFF
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
 	.collectFile(name: "${params.outdir}/Exonerate_output.gff")


/*
 * STEP 6 - Exonerate to Hints
 */
 
process Exonerate2Hints {

	input:
	file exonerate_input from exonerate_for_hints.collectFile()
	
	output:
	file exonerate_hints
	
	"""
	grep -v '#' $exonerate_input | grep -e 'CDS' -e 'exon' -e 'intron' | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > exonerate_hints
	"""
}

exonerate_hints
	.collectFile(name: "${params.outdir}/Exonerate_protein_hints.gff")	


/*
 * STEP 7 - GenomeThreader
 */
 
process RunGenomeThreader {
	
	publishDir "${params.outdir}/genomethreader", mode: 'copy'
		
	output:
	file output_gth
	
	"""
	gth -genomic $Genome -protein $Queries -gff3out -intermediate -o output_gth
	"""
}


/*
 * STEP 8 - GenomeThreader to Hints
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
	.collectFile(name: "${params.outdir}/GenomeThreader_protein_hints.gff")


//RepeatMasker Block

Channel
	.fromPath(Genome)
	.splitFasta(by: params.nrepeats, file: true)
	.ifEmpty { exit 1, "Could not find genome file" }
	.set {fasta_rep}


/*
 * STEP 9 - RepeatMasker
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
 * STEP 10 - RepeatMasker - Collect and Clean1
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
 * STEP 11 - RepeatMasker - Clean2
 */
 
process CleanRepeatMasker {
	
	input:
	file mergedUNCLEAN
	
	output:
	file RepeatMasker_out into RM_clean_out, RM_2_hints
	"""
	grep -v 'with_header' $mergedUNCLEAN | awk 'NF' > RepeatMasker_out
	"""
}

RM_clean_out
 	.collectFile(name: "${params.outdir}/RepeatMasker_output.txt")


/*
 * STEP 12 - RepeatMasker to Hints
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
	.collectFile(name: "${params.outdir}/RepeatMasker_hints.gff")

/*
 * Parse software version numbers
 */
process get_software_versions {

	publishDir "${params.outdir}/report", mode: 'copy'
	
    output:
    file versions

    """
    tblastn -version | head -n 1 > versions
    """
}


workflow.onComplete {
        log.info "========================================="
        log.info "Duration:             $workflow.duration"
        log.info "========================================="
}
