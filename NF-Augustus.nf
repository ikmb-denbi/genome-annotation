#!/usr/bin/env nextflow
/*
========================================================================================
                         NF-Augustus
========================================================================================
 NF-Augustus Analysis Pipeline. Started 2018-08-03.
 #### Homepage / Documentation
 https://git.ikmb.uni-kiel.de/m.torres/NF-hints.git
 #### Authors
 MTorres m.torres <m.torres@ikmb.uni-kiel.de> - https://git.ikmb.uni-kiel.de/m.torres>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     NF-Augustus v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run Augustus.nf --genome 'Genome.fasta' --hints 'Hints/*.gff' --model human --retrain true --newspecies turtle -c config/slurm.config --nthreads 3

    Mandatory arguments:
      --genome                      Genome reference
      --hints						Directory with hints file(s) (e.g. --hints 'Hints/*.gff')
      -c                     	 	Configuration file to use 

    Options:
	  --model						Augustus available profile to use (default = human)
	  --retrain						Whether the available Augustus profile should be retrained with models especific for the species of interest [ true | false (default) ]
	  --newspecies					If '--retrain true', give a name to the new profile to create (default = 'NewSpecies')
	  --naugustus					Chunks (# of scaffolds) to divide Augustus jobs [ default = 1 ]
	  --nthreads					Number of cpus for multi-thread mode [default = 1]
	  --cfg							Config file to use instead of default

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
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

params.model = "human"
params.retrain = false
params.newspecies = "NewSpecies"
params.nthreads = 1
params.name = false
params.outdir = "Augustus_output"
params.cfg = "bin/augustus_default.cfg"
params.naugustus = 1


// Validate inputs
if ( params.genome ){
	Genome = file(params.genome)
    if( !Genome.exists() ) exit 1, "Genome file not found: ${Genome}"
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
                                               

NF-Augustus v${params.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'NF-hints'
summary['Pipeline Version'] = params.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Fasta Ref']    = Genome
summary['Hints']		= params.hints
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


Hints = Channel.fromPath( params.hints )


/*
 * STEP 1 - Concatenate Hints files
 */
process Concatenate {

	publishDir "${params.outdir}", mode: 'copy'

    input:
    file hints2conc from Hints.collect()

    output:
    file 'All_Hints.gff' into all_hints
	
    """
    cat $hints2conc >> All_Hints.gff
    """
}

Channel
		.fromPath(Genome)
		.splitFasta(by: params.naugustus, file: true)
		.set {fasta_aug}
		
/*
 * STEP 2 - Augustus
 */
process runAugustus1 {

	publishDir "${params.outdir}/Augustus_run1/", mode: 'copy'
	
	input:
	file all_hints
	file fasta_aug
	
	
	output:
	file 'Augustus_out' into augustus_out_gff, augustus_out2train
	
	"""
	grep '>' $fasta_aug | perl -ple 's/^>(\S+)\s?.*/\$1/' > scafs
	grep -F -f scafs $all_hints > scafs_hints
	augustus --species=$params.model --UTR=off --alternatives-from-evidence=false --extrinsicCfgFile=$params.cfg --hintsfile=scafs_hints fasta_aug > Augustus_out
	"""
}

augustus_out_gff
	.collectFile( name: "${params.outdir}/Augustus_run1.gff" )


/*
 * STEP 3 - Select Best Models
 */

//process BestModels {

//	publishDir "${params.outdir}/best_models/", mode: 'copy'
	
//	input:
//	file all_models from augustus_out2train.collect()
	
//	output:
	
//	when:
//	params.retrain == true
	
//	"""
	
//	"""
//}





workflow.onComplete {

    log.info "========================================="
    log.info "Duration:             $workflow.duration"
    log.info "========================================="
        
}
