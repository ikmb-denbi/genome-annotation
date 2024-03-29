#!/usr/bin/env nextflow
/*
========================================================================================
                    IKMB - de.NBI | Genome Annotation Pipeline
========================================================================================
 Genome Annotation Pipeline. Started 2018-10-17.
 #### Homepage / Documentation
 https://github.com/ikmb-denbi/genome-annotation/
 #### Authors
 MTorres m.torres <m.torres@ikmb.uni-kiel.de> - https://git.ikmb.uni-kiel.de/m.torres>
 MHoeppner m.hoeppner <m.hoeppner@ikmb.uni-kiel.de> 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  log.info"""
  =================================================================
   IKMB - de.NBI | Genome Annotation Pipeline | v${workflow.manifest.version}
  =================================================================
  Usage:

  The typical command for running the pipeline is as follows:

  nextflow run ikmb-denbi/genome-annotation --genome 'Genome.fasta' --proteins 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' -c nextflow.config

  Mandatory arguments:
  --genome		Genome reference
      
  At least one of:
  --proteins		Proteins from other species
  --ESTs		ESTs or transcriptome
  --reads		Path to RNA-seq data (must be surrounded with quotes)

  Options:
    -profile            Hardware config to use (optional, will default to 'standard')
    Programs to run:
    --trinity		Run transcriptome assembly with Trinity and produce hints from the transcripts [ true (default) | false ]
    --augustus		Run Augustus to predict genes [ true (default) | false ]
    --training		Run de novo model training for Augustus with complete proteins (from Pasa + Transdecoder). Only if RNA-seq data provided. [ true | false ]. You must provide a name for your model with "--model". If you use the name of an existing model, this will be re-trained.
    --pasa 		Run the transcriptome-based gene builder PASA (also required when running --training). [ true | false (default) ]. Requires --ESTs and/or --reads with --trinity. 
 	
    Programs parameters:
    --rm_species	Species database for RepeatMasker [ default = 'mammal' ]
    --rm_lib		Additional repeatmasker library in FASTA format [ default = 'false' ]
    --train_perc	What percentage of complete proteins (from Pasa + Transdecoder) should be used for training. The rest will be used for testing model accuracy [ default = 90 ]]
    --training_models	How many PASA gene models to select for training [ default = 1000 ]
    --model		Species model for Augustus [ default = 'human' ]. If "--training true" and you want to do de novo training, give a NEW name to your species
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --max_intron_size	Maximum length of introns to consider for spliced alignments [ default = 20000 ]
    --evm		Whether to run EvicenceModeler at the end to produce a consensus gene build [true | false (default) ]
    --evm_weights	Custom weights file for EvidenceModeler (overrides the internal default)
  
    Evidence tuning
    --pri_prot		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 5)
    --pri_est		A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 3)
    --pri_rnaseq	A positive number between 1 and 5 - the higher, the more important the hint is for gene calling (default: 4)
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide genome for blastx jobs [ default = 100 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nchunks		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]
    --chunk_size 	Size of sub-regions of the genome on which to run Blastx jobs [ default = 50000 ]

    Other options:
    --singleEnd		Specifies that the input is single end reads [ true | false (default) ]
    --rnaseq_stranded	Whether the RNAseq reads were sequenced using a strand-specific method (dUTP) [ true | false (default) ]
    --outdir		The output directory where the results will be saved [ default = 'output' ]
    -name		Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
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

// -----------------------------
// Validate and set input options
// -----------------------------

OUTDIR = params.outdir

pasa_config = "${baseDir}/assets/pasa/alignAssembly.config"

evm_weights = params.evm_weights ?: "${baseDir}/assets/evm/weights.txt"
EVM_WEIGHTS = file(evm_weights)

uniprot_path = "${baseDir}/assets/Eumetazoa_UniProt_reviewed_evidence.fa"

if (params.pasa) {
	PASA_CONFIG = file(pasa_config)
	if ( !PASA_CONFIG.exists()) exit 1; "Could not find the pasa config file that should be bundled with this pipeline, exiting..."
}

Uniprot = file(uniprot_path)
if ( !Uniprot.exists()) exit 1; "Could not find the Uniprot data that should be bundled with this pipeline, exiting...."

Genome = file(params.genome)
if( !Genome.exists() || !params.genome ) exit 1; "No genome assembly found, please specify with --genome"

if (params.proteins) {
	Proteins = file(params.proteins)
	if( !Proteins.exists() ) exit 1, "Protein file not found: ${Proteins}. Specify with --proteins."
}

if ( params.ESTs ){
	ESTs = file(params.ESTs)
	if( !ESTs.exists() ) exit 1, "ESTs file not found: ${ESTs}. Specify with --ESTs."
}

if (params.rm_lib) {
	RM_LIB = file(params.rm_lib)
	if (!RM_LIB.exists() ) exit 1, "Repeatmask library does not exist (--rm_lib)!"
	if (params.rm_species) {
		println "Provided both a custom repeatmask library (--rm_lib) AND a species/taxonomic group for RM - will only use the library!"
	}
}

// Make it fail if basic requirements are unmet
if (!binding.variables.containsKey("Proteins") && !binding.variables.containsKey("ESTs") && params.reads == false) {
	exit 1, "At least one type of input data must be specified (--proteins, --ESTs, --reads)"
}

if (params.trinity && !params.reads ) {
	exit 1, "Cannot run Trinity de-novo assembly without RNA-seq reads (specify both --reads and --trinity)"
}

// Use a default config file for Augustus if none is provided
if (params.augustus  && !params.augCfg ) {
	AUG_CONF = "$workflow.projectDir/assets/augustus/augustus_default.cfg"
} else if (params.augustus) {
	AUG_CONF = file(params.augCfg)
}

// Check prereqs for training a new model
if (params.training && !params.model) {
        exit 1; "You requested for a new prediction profile to be trained, but did not provide a name for the model (--model)"
} else if (params.training && !params.trinity && !params.ESTs) {
        exit 1; "You requested for a new prediction profile to be trained, but we need transcriptome data for that (--trinity and/or --ESTs)"
} else if (params.training && !params.pasa) {
	println "You requested a model to be trained; this requires Pasa to be enabled. We will do that for your now..."
	params.pasa = true
}

// Check if we can run EVM
if (params.evm && !params.augustus && !params.pasa) {
	exit 1; "Requested to run EvidenceModeler, but we need gene models for that (--augustus and/or --pasa)."
}

//Check if we can run Pasa
if (params.pasa && !params.ESTs && !params.trinity) {
	exit 1; "Requested to run Pasa, but we need transcriptome data for that (--ESTs or --trinity)"
}

// Check prereqs for repeatmasking
if (!params.rm_lib && !params.rm_species) {
	println "No repeat library provided, will model repeats de-novo instead using RepeatModeler."
}

// give this run a name
run_name = ( !params.run_name) ? "${workflow.sessionId}" : "${params.run_name}"

def summary = [:]

summary['Assembly'] = params.genome
summary['ESTs'] = params.ESTs
summary['Proteins'] = params.proteins
summary['RNA-seq'] = params.reads
summary['RM species'] = params.rm_species
summary['RM library'] = params.rm_lib
summary['Augustus model'] = params.model
summary['ModelTraining'] = params.training
// ----------------------
// ----------------------
// Starting the pipeline
// ----------------------
// ----------------------

// --------------------------
// Set various input channels
// --------------------------

Channel.fromPath(Genome)
	.set { GenomeHisat }

// Split the genome for parallel processing
Channel
	.fromPath(Genome)
	.set { genome_for_splitting }

// if proteins are provided
if (params.proteins ) {

	// goes to blasting of proteins
        Channel
        .fromPath(Proteins)
        .set { protein_to_blast_db }

	// create a cdbtools index for the protein file
	Channel
	.fromPath(Proteins)
	.set { index_prots }
} else {
	prot_exonerate_hints = Channel.empty()
	// Protein Exonerate files to EVM
	exonerate_protein_evm = Channel.empty()
}

// if ESTs are provided
if (params.ESTs) {
	// goes to blasting the ESTs
        Channel
                .fromPath(ESTs)
                .set {fasta_ests}

	// create a cdbtools index for the EST file
	Channel
		.fromPath(ESTs)
		.into { ests_index; est_to_pasa }
} else {
	// EST hints to Augustus
	est_minimap_hints = Channel.empty()
	// EST file to Pasa assembly
	est_to_pasa = Channel.empty()
	// EST exonerate files to EVM
	minimap_ests_to_evm = Channel.empty()
}

// if RNAseq reads are provided
if (params.reads) {

	// Make a HiSat index
	Channel
	        .fromPath(Genome)
        	.set { inputMakeHisatdb }

	// Pass reads to trimming
	Channel
		.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
		.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
		.set {read_files_trimming }

	// can use reads without wanting to run a de-novo transcriptome assembly
	if (!params.trinity) {
	        trinity_minimap_hints = Channel.empty()
	 	minimap_trinity_to_evm = Channel.empty()
	        trinity_to_pasa = Channel.empty()
	}

} else {
	// Trinity hints to Augustus
	trinity_minimap_hints = Channel.empty()
	// RNAseq hints to Augustus
	rnaseq_hints = Channel.empty()
	// Trinity assembly to Pasa assembly
	trinity_to_pasa = Channel.empty()
	// Trinity exonerate files to EVM
	minimap_trinity_to_evm = Channel.empty()
}

if (!params.pasa) {
	// Pasa data to EVM
	pasa_genes_to_evm = Channel.empty()
	pasa_align_to_evm = Channel.empty()
}
// Trigger de-novo repeat prediction if no repeats were provided
if (!params.rm_lib  && !params.rm_species ) {
	Channel
		.fromPath(Genome)
		.set { inputRepeatModeler }
} else {
        repeats_fa = Channel.empty()
}

// Provide the path to the augustus config folder
// If it's in a container, use the hard-coded path, otherwise the augustus env variable
if (!workflow.containerEngine) {
	Channel.from(file(System.getenv('AUGUSTUS_CONFIG_PATH')))
		.ifEmpty { exit 1; "Looks like the Augustus config path is not set? This shouldn't happen!" }
        	.set { augustus_config_folder }
} else {
// this is a bit dangerous, need to make sure this is updated when we bump to the next release version
	Channel.from(file("/opt/conda/envs/genome-annotation-1.0/config"))
        	.set { augustus_config_folder }
}

// Header log info
log.info "========================================="
log.info "ESGA Genome Annotation Pipeline v${workflow.manifest.version}"
log.info "Genome assembly: 		${params.genome}"
if (params.rm_lib) {
	log.info "Repeatmasker lib:		${params.rm_lib}"
} else if (params.rm_species) {
	log.info "Repeatmasker species:		${params.rm_species}"
} else {
	log.info "Repeatmasking:			Compute de-novo"
}
log.info "-----------------------------------------"
log.info "Evidences:"
log.info "Proteins:			${params.proteins}"
log.info "ESTs:				${params.ESTs}"
log.info "RNA-seq:			${params.reads}"
if (params.augustus) {
	log.info "Augustus profile		${params.model}"
}
if (params.augustus && AUG_CONF) {
	log.info "Augustus config file		${AUG_CONF}"
}
if (params.training) {
	log.info "Model training:			yes (${params.model})"
}
log.info "-----------------------------------------"
log.info "Parallelization settings"
log.info "Chunk size for assembly:		${params.chunk_size}"
log.info "Chunk size for Blast:			${params.nblast}"
log.info "Chunk size for Exonerate:		${params.nexonerate}"
log.info "Chunk size for RepeatMasker:		${params.nchunks}"
log.info "-----------------------------------------"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="

def check_file_size(fasta) {

        if (file(fasta).isEmpty()) {
                log.info "No repeats were modelled, will use the built-in repeat library that ships with RepeatMasker"
        }
}

// *******************************************
// Split Genome in parts of roughly equal size
// *******************************************

process splitGenome {

	input:
	file(genome_fa) from genome_for_splitting

	output:
	file("*_chunk_*") into fasta_chunks

	script:
	ref_name = genome_fa.getBaseName() + "_chunk_%3.3d"

	"""
		fastasplitn -in $genome_fa -n ${params.nchunks} -t $ref_name
	"""

}

// We split the list of chunks into channel emissions
fasta_chunk_for_rm_lib = fasta_chunks.flatMap()

// ************************************
// Model Repeats if nothing is provided
// ************************************

if (!params.rm_lib && !params.rm_species) {

	process repeatModel {

	        publishDir "${OUTDIR}/repeatmodeler/", mode: 'copy'

		scratch true

		input:
		file(genome_fa) from inputRepeatModeler

		output:
		file(repeats) into (repeats_fa, check_repeats)

		script:	

		repeats = "consensi.fa"	
		"""
			BuildDatabase -name genome_source -engine ncbi $genome_fa
			RepeatModeler -engine ncbi -pa ${task.cpus} -database genome_source
			cp RM_*/consensi.fa . 
		"""
	}

	// Let the user know if the repeat search was unsuccesful
	check_repeats
	        .filter { fasta -> check_file_size(fasta) }
        	.set  { checked_repeats }
}

// ---------------------------
// RUN REPEATMASKER
//----------------------------

// RepeatMasker library needs ot be writable. Need to do this so we can work with locked containers
process repeatLib {

	label 'short_running'

	publishDir "${OUTDIR}/repeatmasker/", mode: 'copy'

	output:
	file("Library") into RMLibPath

	script:
	
	if (params.rm_lib) {
	"""
		cp ${baseDir}/assets/repeatmasker/my_genome.fa .
		cp ${baseDir}/assets/repeatmasker/repeats.fa .
		mkdir -p Library
                cp ${baseDir}/assets/repeatmasker/DfamConsensus.embl Library/
                gunzip -c ${baseDir}/assets/repeatmasker/taxonomy.dat.gz > Library/taxonomy.dat
		export REPEATMASKER_LIB_DIR=\$PWD/Library
		RepeatMasker -lib repeats.fa my_genome.fa > out
	"""	
	} else {

	"""
		mkdir -p Library
		cp ${baseDir}/assets/repeatmasker/DfamConsensus.embl Library/
		gunzip -c ${baseDir}/assets/repeatmasker/taxonomy.dat.gz > Library/taxonomy.dat
	"""
	}
}

// To get the repeat library path combined with each genome chunk, we do this...
// toString() is needed as RM touches the location each time it runs and thus modifies it. 
rm_lib_path = RMLibPath
        .map { it.toString() }
	.combine(fasta_chunk_for_rm_lib)

// generate a soft-masked sequence for each assembly chunk
// if nothing was masked, return the original genome sequence instead and an empty gff file. 
process repeatMask {

	//scratch true

	publishDir "${OUTDIR}/repeatmasker/chunks"

	input: 
	file(repeats) from repeats_fa.collect().ifEmpty('')
	set env(REPEATMASKER_LIB_DIR),file(genome_fa) from rm_lib_path

	output:
	file(genome_rm) into RMFastaChunks
	file(genome_rm) into (genome_to_minimap_chunk,genome_chunk_to_augustus)
	file(rm_gff) into RMGFF
	set file(rm_gff),file(rm_tbl),file(rm_out)

	script:

	def options = ""
	if (params.rm_lib) {
		options = "-lib $params.rm_lib"
	} else if (params.rm_species) {
		options = "-species $params.rm_species"
	} else {
		if (repeats.size() == 0) {
			options = "-species eukaryota"
		} else {
			options = "-lib $repeats"
		}
	}

	genome_rm = "${genome_fa.getName()}.masked"
	rm_gff = "${genome_fa.getName()}.out.gff"
	rm_tbl = "${genome_fa.getName()}.tbl"
	rm_out = "${genome_fa.getName()}.out"
	
	"""
		echo \$REPEATMASKER_LIB_DIR > lib_dir.txt
		RepeatMasker $options -gff -xsmall -q -pa ${task.cpus} $genome_fa

		test -f ${genome_rm} || cp $genome_fa $genome_rm && touch $rm_gff

	"""
}

// Merge the repeat-masked assembly chunks
process repeatMerge {

	label 'short_running'

        publishDir "${OUTDIR}/repeatmasker", mode: 'copy'

	input:
	file(genome_chunks) from RMFastaChunks.collect()

	output:
	file(masked_genome) into (rm_to_blast_db, rm_to_partition)
	set file(masked_genome),file(masked_genome_index) into (RMGenomeIndexProtein, genome_to_trinity_minimap, genome_to_pasa, genome_to_evm, genome_to_evm_merge, RMGenomeMinimapEst, genome_to_minimap_pasa)

	script:
	
	masked_genome = "${Genome.baseName}.rm.fa"
	masked_genome_index = masked_genome + ".fai"

	"""
		cat $genome_chunks >> merged.fa
		fastasort -f merged.fa > $masked_genome
		samtools faidx $masked_genome
		rm merged.fa
	"""	
}

// ---------------------
// PROTEIN DATA PROCESSING
// ---------------------
if (params.proteins) {
	// ----------------------------
	// Protein BLAST against genome
	// ----------------------------

	// Split the genome into smaller chunks for Blastx
	process prepSplitAssembly {

		publishDir "${OUTDIR}/databases/genome/", mode: 'copy'

		label 'short_running'

		input:
		file(genome_fa) from rm_to_blast_db

		output:
		file(genome_chunks) into genome_chunks_blast
		file(genome_agp) into genome_chunks_agp

		script:

		genome_chunks = genome_fa + ".chunk"
		genome_agp = genome_fa + ".agp"

		"""
			chromosome_chunk.pl -fasta_file $genome_fa -size $params.chunk_size
		"""
	}

	genome_chunks_blast_split = genome_chunks_blast.splitFasta(by: params.nblast, file: true)

	// Make a blast database
	process protMakeDB {

		label 'medium_running'

        	publishDir "${OUTDIR}/databases/blast/", mode: 'copy'

	        input:
        	file(protein_fa) from protein_to_blast_db

	        output:
        	file("${dbName}.dmnd") into blast_db_prots

	        script:
        	dbName = protein_fa.getBaseName()
	        """
			diamond makedb --in $protein_fa --db $dbName
        	"""
	}

	// create a cdbtools compatible  index
	// we need this to do very focused exonerate searches later
	process protIndex {

		label 'short_running'

		publishDir "${OUTDIR}/databases/cdbtools/proteins", mode: 'copy'

		input:
		file(fasta) from index_prots

		output:
		set file(fasta),file(protein_index) into ProteinDB

		script:
		protein_index = fasta.getName()+ ".cidx"

		"""
			cdbfasta $fasta 
		"""
	}

	// Blast each genome chunk against the protein database
	// This is used to define targets for exhaustive exonerate alignments
	process protDiamondx {

		publishDir "${OUTDIR}/evidence/proteins/blastx/chunks", mode: 'copy'

		scratch true

		input:
		file(genome_chunk) from genome_chunks_blast_split
		file(blastdb_files) from blast_db_prots.collect()

		output:
		file(protein_blast_report) into ProteinBlastReport

		script:
		db_name = blastdb_files[0].baseName
		chunk_name = genome_chunk.getName().tokenize('.')[-2]
		protein_blast_report = "${genome_chunk.baseName}.blast"
		"""
			diamond blastx --sensitive --threads ${task.cpus} --evalue ${params.blast_evalue} --outfmt ${params.blast_options} --db $db_name --query $genome_chunk --out $protein_blast_report
		"""
	}

	// Parse Protein Blast output for exonerate processing
	process protDiamondToTargets {

		label 'short_running'

	        publishDir "${OUTDIR}/evidence/proteins/tblastn/chunks", mode: 'copy'

		input:
		file(blast_reports) from ProteinBlastReport.collect()
		file(genome_agp) from genome_chunks_agp

		output:
		file(query2target_result_uniq_targets) into query2target_uniq_result_prots
	
		script:
		query_tag = Proteins.baseName
		query2target_result_uniq_targets = "${query_tag}.targets"
	
		"""
			cat $blast_reports > merged.txt
			blast_chunk_to_toplevel.pl --blast merged.txt --agp $genome_agp > merged.translated.txt
			blast2exonerate_targets.pl --infile merged.translated.txt --max_intron_size $params.max_intron_size > $query2target_result_uniq_targets
			rm merged.*.txt
		"""
	}

	// split Blast hits for parallel processing in exonerate
	query2target_uniq_result_prots
		.splitText(by: params.nexonerate, file: true)
		.combine(ProteinDB)
		.set{ query2target_chunk_prots }

	// Run Exonerate on the blast regions
	process protExonerate {

		publishDir "${OUTDIR}/evidence/proteins/exonerate/chunks", mode: 'copy'

		//scratch true

		input:
		set file(hits_chunk),file(protein_db),file(protein_db_index) from query2target_chunk_prots
		set file(genome),file(genome_faidx) from RMGenomeIndexProtein
	
		output:
		file(exonerate_chunk) into (exonerate_result_prots, exonerate_protein_chunk_evm)
		file(commands)

		script:
		query_tag = protein_db.baseName
		chunk_name = hits_chunk.getName().tokenize('.')[-2]
		commands = "commands." + chunk_name + ".txt"
		exonerate_chunk = "${hits_chunk.baseName}.${query_tag}.exonerate.out"
		
		// get the protein fasta sequences, produce the exonerate command and genomic target interval fasta, run the whole thing,
		// merge it all down to one file and translate back to genomic coordinates
		// remove all the untracked intermediate files

		"""
			extractMatchTargetsFromIndex.pl --matches $hits_chunk --db $protein_db_index
			exonerate_from_blast_hits.pl --matches $hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --query_index $protein_db_index --analysis protein2genome --outfile $commands
			parallel -j ${task.cpus} < $commands
			echo "Finished exonerate run" >> log.txt
			echo "# from ${chunk_name}" > merged.${chunk_name}.exonerate.out ;
			cat *.exonerate.align | grep -v '#' | grep 'exonerate:protein2genome:local' >> merged.${chunk_name}.exonerate.out 2>/dev/null
			echo "Finished merging exonerate alignments" >> log.txt
			exonerate_offset2genomic.pl --infile merged.${chunk_name}.exonerate.out --outfile $exonerate_chunk
			echo "Finished translating genomic coordinates" >> log.txt
			test -f $exonerate_chunk || cat "#" > $exonerate_chunk
		"""
	}

	exonerate_protein_evm = exonerate_protein_chunk_evm.collectFile()

	// merge the exonerate hits and create the hints
	process protExonerateToHints {

		label 'medium_running'

		publishDir "${OUTDIR}/evidence/proteins/exonerate/", mode: 'copy'

		input:
		file(chunks) from exonerate_result_prots.collect()

		output:
		file(exonerate_gff) into prot_exonerate_hints

		script:
		query_tag = Proteins.baseName
		exonerate_gff = "proteins.exonerate.${query_tag}.hints.gff"
		"""
			cat $chunks > all_chunks.out
			exonerate2gff.pl --infile all_chunks.out --pri ${params.pri_prot} --source protein --outfile $exonerate_gff
		"""
	}

} // close protein loop

// --------------------------
// -------------------
// EST DATA PROCESSING
//-------------------
// --------------------------

if (params.ESTs) {

	/*
	 * EST alignment
	*/

	// Align all ESTs against the masked genome
	process estMinimap {

		publishDir "${OUTDIR}/evidence/EST/minimap", mode: 'copy'

		scratch true

		input:
		file(est_chunks) from fasta_ests
		set file(genome_rm),file(genome_index) from RMGenomeMinimapEst

		output:
		file(minimap_gff) into (minimap_est_gff, minimap_ests_to_evm)

		script:
		minimap_gff = "ESTs.minimap.gff"
		minimap_bam = "ESTS.minimap.bam"
		"""
			minimap2 -t ${task.cpus} -ax splice:hq -c -G ${params.max_intron_size}  $genome_rm $est_chunks | samtools sort -O BAM -o $minimap_bam
			minimap2_bam2gff.pl $minimap_bam > $minimap_gff
		"""	

	}

	// Combine exonerate hits and generate hints
	process estMinimapToHints {

		label 'short_running'

		publishDir "${OUTDIR}/evidence/EST/minimap/", mode: 'copy'
	
		input:
		file(minimap_gff) from minimap_est_gff
	
		output:
		file(minimap_hints) into est_minimap_hints
	
		script:
		minimap_hints = "ESTs.minimap.hints.gff"
			
		"""
			minimap2hints.pl --source est2genome --pri ${params.pri_est} --infile $minimap_gff --outfile $minimap_hints
		"""
	}

} // close EST loop

// ++++++++++++++++++
// RNA-seq PROCESSING
// ++++++++++++++++++
if (params.reads) {

	// trim reads
	process rseqTrim {

		publishDir "${OUTDIR}/evidence/rnaseq/fastp", mode: 'copy'

		scratch true 

		input:
		set val(name), file(reads) from read_files_trimming

		output:
		file("*_trimmed.fastq.gz") into trimmed_reads
		set file(json),file(html) into trimmed_reads_qc

		script:
		prefix = reads[0].toString().split("_R1")[0]
		json = prefix + ".fastp.json"
		html = prefix + ".fastp.html"

		if (params.singleEnd) {
			left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
			"""
                       		fastp -i ${reads[0]} --out1 ${left} -w ${task.cpus} -j $json -h $html
	                """
		} else {
			left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
			right = file(reads[1]).getBaseName() + "_trimmed.fastq.gz"
			"""
				fastp --detect_adapter_for_pe --in1 ${reads[0]} --in2 ${reads[1]} --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html
			"""
		}
	}

	// Generate an alignment index from the genome sequence
	process rseqMakeDB {

		label 'long_running'

		publishDir "${OUTDIR}/databases/HisatDB", mode: 'copy'

		input:
		file(genome) from inputMakeHisatdb
	
		output:
		file "${dbName}.*.ht2" into hs2_indices
	
		script:
		dbName = genome.baseName
		dbName_1 = dbName + ".1.ht2"
		
		prefix = dbName
		"""
		hisat2-build $genome $dbName -p ${task.cpus}
		"""
	}

	/*
	 * STEP RNAseq.3 - Hisat2
	 */

	process rseqMap {

		publishDir "${OUTDIR}/evidence/rnaseq/Hisat2/libraries", mode: 'copy'
	
		scratch true

		input:
		file reads from trimmed_reads
		file hs2_indices from hs2_indices.collect()	
	
		output:
		file "*accepted_hits.bam" into accepted_hits2merge
	
		script:
		indexBase = hs2_indices[0].toString() - ~/.\d.ht2/
		ReadsBase = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

		prefix = ReadsBase

		if (params.singleEnd) {
			"""
			hisat2 -x $indexBase -U $reads -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@ 4 - > ${prefix}_accepted_hits.bam
			"""
		} else {
			"""
			hisat2 -x $indexBase -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} | samtools view -bS - | samtools sort -m 2G -@4 - > ${prefix}_accepted_hits.bam
			"""
   		}
	}

	// Combine all BAM files for hint generation
	process rseqMergeBams {

		publishDir "${OUTDIR}/evidence/rnaseq/Hisat2", mode: 'copy'

		scratch true 

		input:
		file hisat_bams from accepted_hits2merge.collect()

		output:
		file(bam) into  (Hisat2Hints, bam2trinity)

		script:
		bam = "hisat2.merged.bam"
		avail_ram_per_core = (task.memory/task.cpus).toGiga()-1
	
		"""
			samtools merge - $hisat_bams | samtools sort -@ ${task.cpus} -m${avail_ram_per_core}G - > $bam
		"""
	}

	/*
	 * STEP RNAseq.4 - Hisat2 into Hints
	 */	
	process rseqHints {

		publishDir "${OUTDIR}/evidence/rnaseq/hints", mode: 'copy'

		input:
		file(bam) from Hisat2Hints
	
		output:
		file(hisat_hints) into rnaseq_hints
	
		script:
		hisat_hints = "RNAseq.hisat.hints.gff"

		"""
			bam2hints --intronsonly 0 -p ${params.pri_rnaseq} -s 'E' --in=$bam --out=$hisat_hints
		"""
	}

	// ----------------------------
	// run trinity de-novo assembly
	// ----------------------------

	if (params.trinity) {

		process rseqTrinity {
	
			publishDir "${OUTDIR}/evidence/rnaseq/trinity", mode: 'copy'

			scratch true 
	
			input:
			file(hisat_bam) from bam2trinity

			output:
			file "transcriptome_trinity/Trinity-GG.fasta" into (trinity_transcripts, trinity_to_minimap , trinity_to_pasa)
	
			script:

			trinity_option = ( params.rnaseq_stranded == true ) ? "--SS_lib_type RF" : ""

			"""
				Trinity --genome_guided_bam $hisat_bam \
				--genome_guided_max_intron ${params.max_intron_size} \
				--CPU ${task.cpus} \
				--max_memory ${task.memory.toGiga()-1}G \
				--output transcriptome_trinity \
				$trinity_option
			"""
		}

		process rseqMinimapTrinity {

			publishDir "${OUTDIR}/evidence/rnaseq/trinity", mode: 'copy'

			input:
			file(trinity_fasta) from trinity_to_minimap
			set file(genome_rm),file(genome_index) from genome_to_trinity_minimap

			output:
			file(trinity_gff) into (minimap_trinity_gff,minimap_trinity_to_evm)
			

			script:
			trinity_gff = "trinity.minimap.gff"

			"""
				minimap2 -t ${task.cpus} -ax splice:hq -c $genome_rm $trinity_fasta | samtools sort -O BAM -o minimap.bam
				minimap2_bam2gff.pl minimap.bam > $trinity_gff
			"""
		}

		process rseqTrinityToHints {

			label 'short_running'

			publishDir "${OUTDIR}/evidence/EST/minimap/", mode: 'copy'

	                input:
        	        file(trinity_gff) from minimap_trinity_gff

                	output:
	                file(minimap_hints) into trinity_minimap_hints

        	        script:
                	minimap_hints = "trinity.minimap.hints.gff"

	                """
        	                minimap2hints.pl --source est2genome --pri ${params.pri_est} --infile $trinity_gff --outfile $minimap_hints
                	"""
		}


	} // Close Trinity loop
	
} // Close RNAseq loop

/*
* RUN PASA WITH TRANSCRIPTS
*/

if (params.pasa) {

	// use trinity transcripts, ESTs or both
	if (params.trinity != false || params.ESTs != false ) {

		// Clean transcripts
		// This is a place holder until we figure out how to make Seqclean work within Nextflow
		process transSeqclean {

			label 'short_running'

			publishDir "${OUTDIR}/evidence/rnaseq/pasa/seqclean/", mode: 'copy'

			input:
			file(trinity) from trinity_to_pasa.ifEmpty('')
			file(ests) from est_to_pasa.ifEmpty('')

			output:
			set file(transcripts_clean),file(transcripts) into (seqclean_to_pasa,seqclean_to_minimap)

			script:
			transcripts = "transcripts.fa"
			transcripts_clean = transcripts + ".clean"

			file_list = ""
			if (params.ESTs) {
				file_list += ests
			}
			if (params.trinity) {
				file_list += " ${trinity}"
			}

			"""
				cat $file_list | grep -v false >> $transcripts
				cp $transcripts $transcripts_clean

			"""		
		}

		// run minimap for fast transcript mapping
		process transPasaMinimap {

			publishDir "${OUTDIR}/evidence/transcripts/minimap", mode: 'copy'
		
			input:
			set file(transcripts_clean),file(transcripts) from seqclean_to_minimap
			set file(genome),file(genome_index) from genome_to_minimap_pasa
			output:
			set file(transcripts_clean),file(minimap_gff) optional true into minimap_to_pasa

			script:
			minimap_gff = "minimap.transcripts.gff"	
			minimap_bam = "minimap.transcripts.bam"

			"""
				minimap2 -t ${task.cpus} -ax splice:hq -c $genome $transcripts_clean  | samtools sort -O BAM -o $minimap_bam
				minimap2_bam2gff.pl $minimap_bam > out.gff
				if [ -s out.gff ] ; then 
					mv out.gff $minimap_gff
				fi 
					
			"""
	
		}
		
		// We parallelize PASA by filtering the minimap alignments per genome chunk
		process transPasaMinimapSplit {
	
			label 'short_running'

			input:
			file(genome_chunk) from genome_to_minimap_chunk
			set file(transcripts),file(minimap_gff) from minimap_to_pasa.collect()
			
			output:
			set file(genome_chunk),file(transcripts_minimap),file(minimap_chunk) into minimap_chunk_to_pasa

			script:
			minimap_chunk = genome_chunk.getBaseName() + ".minimap.gff"
			transcripts_minimap = genome_chunk.getBaseName() + ".transcripts.fasta"
			genome_chunk_index = genome_chunk + ".fai"
			// filter the gff file to only contain entries for our scaffolds of interest
			// then make a list of all transcript ids and extract them from the full transcript fasta
			"""
				samtools faidx $genome_chunk
				minimap_filter_gff_by_genome_index.pl --index $genome_chunk_index --gff $minimap_gff --outfile  $minimap_chunk
				minimap_gff_to_accs.pl --gff $minimap_chunk | sort -u > list.txt
				faSomeRecords $transcripts list.txt $transcripts_minimap
			"""

		}

		// Run the PASA pipeline
		process transPasa {
			
			scratch true
			
			publishDir "${OUTDIR}/annotation/pasa/models", mode: 'copy'

			input:
			set file(genome_rm),file(transcripts_minimap),file(minimap_chunk_gff) from minimap_chunk_to_pasa

			output:
			set file(pasa_assemblies_fasta),file(pasa_assemblies_gff) into PasaResults
	
			script:
			trunk = genome_rm.getBaseName()
			pasa_assemblies_fasta = "pasa_DB_${trunk}.sqlite.assemblies.fasta"
			pasa_assemblies_gff = "pasa_DB_${trunk}.sqlite.pasa_assemblies.gff3"
			
			// optional MySQL support
			// create a config file with credentials, add config file to pasa execute
			// and somehow check that the DB doesn't already exists
			mysql_create_options = ""
			mysql_config_option = ""
			mysql_db_name = ""
			if (params.pasa_mysql_user) {
				mysql_options = "make_pasa_mysql_config.pl --infile \$PASAHOME/pasa_conf/conf.txt --outfile pasa_mysql_conf.txt --user ${params.pasa_mysql_user} --pass ${params.pasa_mysql_pass} --host ${params.pasa_mysql_host} --port ${params.pasa_mysql_port}"	
				mysql_config_option = "-C pasa_mysql_conf.txt"
				mysql_db_name = "--mysql $run_name"
			}

			// The pasa sqlite file must have a fully qualified path, the script is a workaround as this seems difficult to do inside 
			// the script statement

			"""
				make_pasa_config.pl --infile $pasa_config --trunk $trunk --outfile pasa_DB.config $mysql_db_name
				$mysql_create_options
				\$PASAHOME/Launch_PASA_pipeline.pl \
					-c pasa_DB.config -C -R \
					-t $transcripts_minimap \
					-I $params.max_intron_size \
					-g $genome_rm \
					--IMPORT_CUSTOM_ALIGNMENTS_GFF3 $minimap_chunk_gff \
					--CPU ${task.cpus} \
					$mysql_config_option
			"""	
		}

		// Extract gene models from PASA database
		// All chunks are merged using a perl script into the base name pasa_db_merged
		// this is not...ideal. 
		process transPasaToModels {

			scratch true

			label 'long_running'

	                publishDir "${OUTDIR}/annotation/pasa/", mode: 'copy'

			input:
                        file(pasa_assemblies) from PasaResults.collect()

			output:
			file(pasa_transdecoder_fasta) 
			file (pasa_transdecoder_gff) into (pasa_to_training, pasa_genes_to_evm)
			file(merged_gff) into pasa_align_to_evm

			script:
			base_name = "pasa_db_merged"
			merged_fasta = base_name + ".assemblies.fasta"
			merged_gff = base_name + ".pasa_assemblies.gff"
			pasa_transdecoder_fasta = merged_fasta + ".transdecoder.pep"
			pasa_transdecoder_gff = merged_fasta + ".transdecoder.genome.gff3"

			script:

			"""
				pasa_merge_chunks.pl --base $base_name

				\$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi \
				--pasa_transcripts_fasta $merged_fasta \
				--pasa_transcripts_gff3 $merged_gff \
                        	
			"""
		}
		
		if (params.training) {
			// Extract full length models for training
			process transModelsToTraining {

				label 'short_running'

                		publishDir "${OUTDIR}/annotation/pasa/training", mode: 'copy'

				input:
				file(pasa_transdecoder_gff) from pasa_to_training
			
				output:
				file(training_gff) into models2train

				script:
				training_gff = "transdec.complete.gff3"
	
				"""
					pasa_select_training_models.pl --nmodels $params.training_models --infile $pasa_transdecoder_gff >> $training_gff
				"""
			}

			// If we have to modify the AUGUSTUS config folder, we must copy it first
		        process prepAugustusConfig {

				label 'short_running'

				publishDir "${OUTDIR}/augustus/", mode: 'copy'
	
        		        input:
                		val(config_folder) from augustus_config_folder

	                	output:
	        	        file(copied_folder) into (acf_training,acf_training_path)

        	        	script:
		                copied_folder = "config_folder"
	
        		        """
                		        cp -R $config_folder $copied_folder
	                	"""
        		}

			// Run one of two training routines for model training
			process trainAugustus {

				label 'extra_long_running'

				scratch true
	
				publishDir "${OUTDIR}/augustus/training/", mode: 'copy'

				input:
				file(complete_models) from models2train
				env AUGUSTUS_CONFIG_PATH from acf_training
				val(acf_folder) from acf_training_path

				output:
				file(training_stats)
				val(acf_folder) into acf_prediction

				script:
				complete_gb = "complete_peptides.raw.gb"	
				train_gb = "complete_peptides.raw.gb.train"
				test_gb = "complete_peptides.raw.gb.test"
				training_stats = "training_accuracy.out"

				// If the model already exists, do not run new_species.pl
				//model_path = "${acf_training_path}/species/${params.model}"
				model_file = file("${acf_folder}/species/${params.model}")
			
				options = ""
				if (!model_file.exists()) {
					options = "new_species.pl --species=${params.model}"
				}

				"""
					echo ${acf_folder.toString()} >> training.txt
					gff2gbSmallDNA.pl $complete_models $Genome 1000 $complete_gb
					split_training.pl --infile $complete_gb --percent 90
					$options
					etraining --species=$params.model --stopCodonExcludedFromCDS=false $train_gb
					optimize_augustus.pl --species=$params.model $train_gb --cpus=${task.cpus} --UTR=off 
					augustus --stopCodonExcludedFromCDS=false --species=$params.model $test_gb | tee $training_stats
				"""
			} 
		} else {
			// Make sure we set the built in augustus config path if no training was done
		        acf_prediction = augustus_config_folder
		}// close training loop
	} // close pasa loop
} else {
	// We have to make the AUGUSTUS_CONFIG_PATH a mutable object, so we have to carry through the location of the modifiable 
	// copy of the original config folder. 
        acf_prediction = augustus_config_folder
	pasa_output = Channel.empty()
}

// get all available hints and merge into one file
process prepMergeHints {

	label 'short_running'

	publishDir "${OUTDIR}/evidence/hints", mode: 'copy'

	input:
	file(protein_exonerate_hint) from prot_exonerate_hints.ifEmpty('')
	file(rnaseq_hint) from rnaseq_hints.ifEmpty('')
        file(est_minimap_hint) from est_minimap_hints.ifEmpty('')
        file(trinity_minimap_hint) from trinity_minimap_hints.ifEmpty('')

	output:
	file(merged_hints) into (mergedHints,mergedHintsSort)

	script:
	def file_list = ""
	file_list += !protein_exonerate_hint.empty() ? " ${protein_exonerate_hint}" : ""
	file_list += !rnaseq_hint.empty() ? " ${rnaseq_hint}" : ""
	file_list += !est_minimap_hint.empty() ? " ${est_minimap_hint}" : ""
	file_list += !trinity_minimap_hint.empty() ? " ${trinity_minimap_hint}" : ""

	merged_hints = "merged.hints.gff"
	
	"""
		cat $file_list | grep -v "false" >> $merged_hints
	"""
}

// ****************************
// Split hints by strand and calculate windows for downstream annotation
// ****************************
process prepHintsToBed {

	publishDir "${OUTDIR}/annotation/augustus/targets"

	label 'short_running'

	input:
	file(hints) from mergedHintsSort

	output:
	file(bed) into HintRegions

	script:

	bed = "regions.bed"

	"""
                grep -v "#" $hints | grep -v "false" | sort -k1,1 -k4,4n -k5,5n -t\$'\t' > hints.sorted

		gff_by_strand.pl --infile hints.sorted
		gff2clusters.pl --infile hints.plus.gff --max_intron $params.max_intron_size > $bed
                gff2clusters.pl --infile hints.minus.gff --max_intron $params.max_intron_size >> $bed

	"""
}

// execute Augustus
/*
 * STEP Augustus.1 - Genome Annotation
 */
// Run against each repeatmasked chunk of the assembly

augustus_input_chunk = acf_prediction
	.map { it.toString() }
	.combine(genome_chunk_to_augustus)

process predAugustus {

	publishDir "${OUTDIR}/annotation/augustus/chunks"

	//scratch true 

        when:
        params.augustus

	input:
	set env(AUGUSTUS_CONFIG_PATH),file(genome_chunk) from augustus_input_chunk
	file(regions) from HintRegions.collect()
	file(hints) from mergedHints.collect()

	output:
	file(augustus_result) into augustus_out_gff
	file(command_file)

	script:
	chunk_name = genome_chunk.getName().tokenize(".")[-2]
	augustus_result = "augustus.${chunk_name}.out.gff"
	genome_fai = genome_chunk.getName() + ".fai"
	command_file = "commands." + chunk_name + ".txt"

	"""
		samtools faidx $genome_chunk
		fastaexplode -f $genome_chunk -d . 
		augustus_from_regions.pl --genome_fai $genome_fai --model $params.model --utr ${params.utr} --isof false --aug_conf $AUG_CONF --hints $hints --bed $regions > $command_file
		parallel -j ${task.cpus} < $command_file
		touch dummy.augustus.gff
		cat *augustus.gff > $augustus_result
	"""
}

// Merge all the chunk GFF files into one file
process prepMergeAugustus {

	label 'short_running'

	publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'
	
	input:
	file(augustus_gffs) from augustus_out_gff.collect()

	output:
	file(augustus_merged_gff) into (augustus_2prots, augustus_gff2functions, augustus_to_evm)

	script:
	augustus_merged_gff = "augustus.merged.out.gff"
	
	"""	
		cat $augustus_gffs >> merged.gff
		create_gff_ids.pl --gff merged.gff > $augustus_merged_gff
	"""
}

// Dump out proteins for subsequent functional annotation
process prepAugustusToProtein {

	label 'short_running'

        publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'

	input:
	file augustus2parse from augustus_2prots
	
	output:
	file(augustus_prot_fa) into (augustus_prots2annie, augustus_prots2functions)
	
	script:
	augustus_prot_fa = "augustus.proteins.fa"

	"""
		getAnnoFasta.pl $augustus2parse
		cat *.aa > $augustus_prot_fa	
	"""
}

/*
 * Run EvidenceModeler
*/

if (params.evm) {

	process predEvmPartition {

		label 'long_running'

		// CANNOT USE SCRATCH FOR THIS!
		//scratch true

		//publishDir "${OUTDIR}/annotation/evm/jobs", mode: 'copy'

		input:
		file(augustus_gff) from augustus_to_evm.ifEmpty(false)
		file(est) from minimap_ests_to_evm.ifEmpty(false)
		file(pasa_align) from pasa_align_to_evm.ifEmpty(false)
		file(trinity) from minimap_trinity_to_evm.ifEmpty(false)
		file(proteins) from exonerate_protein_evm.ifEmpty(false)
		file(pasa_genes) from pasa_genes_to_evm.ifEmpty(false)
		set file(genome_rm),file(genome_index) from genome_to_evm

		output:
		file(evm_commands) into inputToEvm
		file(partitions) into EvmPartition

		script:

		partitions = "partitions_list.out"
		evm_commands = "commands.list"
		transcripts = "transcripts.merged.gff"
                gene_models = "gene_models.gff"

		protein_options = ""
		transcript_options = ""
		if (proteins == "false" || proteins == false || proteins =~ /input.*/) {
		} else {
			protein_options = "--protein_alignments $proteins "
		}
		if ( est || trinity || pasa_align ) {
			transcript_options = "--transcript_alignments $transcripts "
		}

		"""
			cat $est $trinity $pasa_align | grep -v false >> $transcripts
			cat $augustus_gff $pasa_genes | grep -v false >> $gene_models

			\$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome $genome_rm \
				--gene_predictions $gene_models \
				--segmentSize 2000000 --overlapSize 200000 --partition_listing $partitions \
				$protein_options $transcript_options
				
			\$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome $genome_rm \
				--weights $EVM_WEIGHTS \
				--gene_predictions $gene_models \
				--output_file_name evm.out \
				--partitions $partitions > $evm_commands
		"""	
			
	}

	// run n commands per chunk
	evm_command_chunks = inputToEvm.splitText(by: params.nevm, file: true)

	// The outputs doesn't do nothing here, EVM combines the chunks based on the original partition file
	process predEvm {

		scratch true

		label 'long_running'

                //publishDir "${OUTDIR}/annotation/evm/chunks", mode: 'copy'

		input:
		file(evm_chunk) from evm_command_chunks

		output:
		file(log_file) into EvmOut

		script:
		log_file = evm_chunk.getBaseName() + ".log"
		"""
			\$EVM_HOME/EvmUtils/execute_EVM_commands.pl $evm_chunk | tee $log_file
		"""
		
	}
	
	// Merge all the separate partition outputs into a final gff file
	process predEvmMerge {
	
		label 'long_running'

		publishDir "${OUTDIR}/annotation/evm", mode: 'copy'

		input:

		file(evm_logs) from EvmOut.collect()
		file(partitions) from EvmPartition
		set file(genome_rm),file(genome_index) from genome_to_evm_merge

		output:
		file(partitions) into EvmResult
		file(done)

		script:
		evm_final = "evm.out"
		done = "done.txt"
		"""
			\$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions $partitions --output_file_name evm.out
			\$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions $partitions --output evm.out --genome $genome_rm
			touch $done
		"""

	}

	// We merge the partial gffs from the partitions with a perl script, 
	// since the output folders are not transferred between processes
	process predEvmToGff {

		label 'medium_running'

		publishDir "${OUTDIR}/annotation/evm", mode: 'copy'

		input:
		file(partitions) from EvmResult

		output:
		file(evm_gff) into EvmGFF

		script:
		evm_gff = "annotations.evm.gff"

		"""
			merge_evm_gff.pl --partitions $partitions --gff $evm_gff
		"""

	}

} // end evm loop


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  summary["Output"] = workflow.launchDir + "/" + OUTDIR
  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "Annotation run finished ($run_name)."

  if (params.email) {

  	def mqc_report = null

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}
