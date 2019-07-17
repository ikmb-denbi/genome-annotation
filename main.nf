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

  nextflow run main.nf --genome 'Genome.fasta' --prots 'Proteins.fasta' --reads 'data/*_R{1,2}.fastq' -c nextflow.config

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
 	
    Programs parameters:
    --rm_species		Species database for RepeatMasker [ default = 'mammal' ]
    --rm_lib		Additional repeatmasker library in FASTA format [ default = 'false' ]
    --model		Species model for Augustus [ default = 'human' ]
    --UTR		Allow Augustus to predict UTRs (results are not optimal and takes much longer) [ 'on' | 'off' (default) ]
    --iso		Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer) [ 'true' | 'false' (default) ]
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --max_intron_size	Maximum length of introns to consider for spliced alignments [ default = 20000 ]
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide Blast jobs [ default = 500 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nrepeats		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]

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

Genome = file(params.genome)
if( !Genome.exists() || !params.genome ) exit 1; "No genome assembly found, please specify with --genome"

if (params.proteins) {
	Proteins = file(params.proteins)
	if( !Proteins.exists() ) exit 1, "Protein file not found: ${Proteins}. Specify with --proteins."
}

if ( params.ESTs ){
	ESTs = file(params.ESTs)
	if( !ESTs.exists() ) exit 1, "ESTs file not found: ${ESTs}. Specify with --ESTs."
	println "Will run Exonerate on EST/Transcriptome file"
}

if (params.reads){
	println "Found reads and will run Hisat2 on RNA-seq data"
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

if (params.trinity == true && params.reads == false) {
	exit 1, "Cannot run Trinity de-novo assembly without RNA-seq reads (specify both --reads and --trinity)"
}

// Use a default config file for Augustus if none is provided
if (params.augustus != false && params.augCfg == false ) {
	AUG_CONF = "$workflow.projectDir/bin/augustus_default.cfg"
	println "Using Augustus config bundled with this pipeline..."
} else if (params.augustus != false) {
	AUG_CONF = params.augCfg
}

if (params.rm_lib == false && params.rm_species == false) {
	println "No repeat library provided, will model repeats de-novo instead using RepeatModeler."
}

// give this run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

def summary = [:]

summary['Assembly'] = params.genome
summary['ESTs'] = params.ESTs
summary['Proteins'] = params.proteins
summary['RNA-seq'] = params.reads
summary['RM species'] = params.rm_species
summary['RM library'] = params.rm_lib
summary['Augustus model'] = params.model

// ----------------------
// ----------------------
// Starting the pipeline
// ----------------------
// ----------------------
// Logic:
// - repeatmask the genome
// -- make masked blast database
// --- align proteins with blast
// ---- align proteins with exonerate
// --- align ESTs with Blast
// -- Align RNAseq reads with HiSat2
// -- Assembly Trinity transcripts (genome-guided)
// -- Generate hints from all available sources
// -- Run AUGUSTUS prediction on each chunk of the masked genome using all available hints

// --------------------------
// Set various input channels
// --------------------------

Channel.fromPath(Genome)
	.set { GenomeHisat }

// Split the genome for parallel processing
Channel
	.fromPath(Genome)
	.splitFasta(by: params.nrepeats, file: true)
	.set { FastaRM }

// if proteins are provided
if (params.proteins) {

	// goes to blasting of proteins
        Channel
        .fromPath(Proteins)
        .splitFasta(by: params.nblast, file: true)
        .set { fasta_prots }

	// create a cdbtools index for the protein file
	Channel
	.fromPath(Proteins)
	.set { index_prots }
} else {
	prot_exonerate_hints = Channel.from(false)
}

// if ESTs are provided
if (params.ESTs) {

	// goes to blasting the ESTs
        Channel
                .fromPath(ESTs)
                .splitFasta(by: params.nblast, file: true)
                .set {fasta_ests}

	// create a cdbtools index for the EST file
	Channel
		.fromPath(ESTs)
		.set { ests_index }
} else {
	est_exonerate_hints = Channel.from(false)
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
	if (params.trinity == false) {
	        trinity_exonerate_hints = Channel.from(false)
	}

} else {
	trinity_exonerate_hints = Channel.from(false)
	rnaseq_hints = Channel.from(false)
}

// Trigger de-novo repeat prediction of no repeats were provided
if (params.rm_lib == false && params.rm_species == false) {
	Channel
		.fromPath(Genome)
		.set { inputRepeatModeler }
} else {
        repeats_fa = Channel.from(false)
}

// Header log info
log.info "========================================="
log.info "IKMB Genome Annotation Pipeline v${workflow.manifest.version}"
log.info "Genome assembly: 		${params.genome}"
if (params.rm_lib) {
	log.info "Repeatmasker lib:		${params.rm_lib}"
} else if (params.rm_species) {
	log.info "Repeatmasker species:		${params.rm_species}"
} else {
	log.info "Repeamasking:			Compute de-novo"
}
log.info "-----------------------------------------"
log.info "Evidences:"
log.info "Proteins:			${params.proteins}"
log.info "ESTs:				${params.ESTs}"
log.info "RNA-seq:			${params.reads}"
if (params.augustus) {
	log.info "Augustus profile		${params.model}"
}
if (params.augCfg) {
	log.info "Augustus config file		${AUG_CONF}"
}
log.info "-----------------------------------------"
log.info "Parallelization settings"
log.info "Chunk size for Blast:			${params.nblast}"
log.info "Chunk size for Exonerate:		${params.nexonerate}"
log.info "Chunk size for RepeatMasker:		${params.nrepeats}"
log.info "-----------------------------------------"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="

// Model Repeats if nothing is provided
def check_file_size(fasta) {

        if (file(fasta).isEmpty()) {
                log.info "No repeats were modelled, will use the built-in repeat library that ships with RepeatMasker"
        }
}

if (params.rm_lib == false && params.rm_species == false ) {

	process runRepeatModeler {

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

// RepeatMasker library needs ot be writable. Need to do this so we can work with locked Singularity containers
process createRMLib {

	tag "ALL"
	publishDir "${OUTDIR}/repeatmasker/", mode: 'copy'

	output:
	file("Library") into RMLibPath

	script:

	"""
		mkdir -p Library
		cp ${baseDir}/assets/repeatmasker/DfamConsensus.embl Library/ 
		gunzip -c ${baseDir}/assets/repeatmasker/taxonomy.dat.gz > Library/taxonomy.dat
	"""

}

// ---------------------------
// RUN REPEATMASKER
//----------------------------

// generate a soft-masked sequence for each assembly chunk
// toString is needed because RM modifies the library each time it touches it.
// if nothing was masked, return the original genome sequence instead and an empty gff file. 
process runRepeatMasker {

	publishDir "${OUTDIR}/repeatmasker/chunks"

	input: 
	file(genome_fa) from FastaRM
	file(repeats) from repeats_fa
	env(REPEATMASKER_LIB_DIR) from RMLibPath.map { it.toString() } 

	output:
	file(genome_rm) into RMFastaChunks
	file(rm_gff) into RMGFF

	script:

	def options = ""
	if (params.rm_lib) {
		options = "-lib $params.rm_lib"
	} else if (params.rm_species) {
		options = "-species $params.rm_species"
	} else {
		if (repeats.size() == 0) {
		} else {
			options = "-lib $repeats"
		}
	}

	genome_rm = "${genome_fa.getName()}.masked"
	rm_gff = "${genome_fa.getName()}.out.gff"
	
	"""
		RepeatMasker $options -gff -xsmall -q -pa ${task.cpus} $genome_fa

		test -f ${genome_rm} || cp $genome_fa $genome_rm && touch $rm_gff

	"""
}

// Merge the repeat-masked assembly chunks
process runMergeRMGenome {

	tag "ALL"
        publishDir "${OUTDIR}/repeatmasker", mode: 'copy'

	input:
	file(genome_chunks) from RMFastaChunks.collect()

	output:
	file(masked_genome) into (RMtoBlastDB, RMtoSplit,RMtoPartition)
	set file(masked_genome),file(masked_genome_index) into RMGenomeIndexProtein, RMGenomeIndexEST, RMGenomeIndexTrinity

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
GenomeChunksAugustus = RMtoPartition
	.splitFasta(by: params.nrepeats, file: true)

// Turn genome into a masked blast database
// Generates a dust mask from softmasked genome sequence
process runMakeBlastDB {
	
	tag "ALL"
	publishDir "${OUTDIR}/databases/blast/", mode: 'copy'

	input:
	file(genome_fa) from RMtoBlastDB

	output:
	file("${dbName}*.n*") into (blast_db_prots, blast_db_ests, blast_db_trinity)

	script:
	dbName = genome_fa.baseName
	
	"""
		makeblastdb -in $genome_fa -parse_seqids -dbtype nucl -out $dbName 
	"""
}

// ---------------------
// PROTEIN DATA PROCESSING
// ---------------------
if (params.proteins != false ) {
	// ----------------------------
	// Protein BLAST against genome
	// ----------------------------

	// create a cdbtools compatible  index
	// we need this to do very focused exonerate searches later
	process runIndexProteinDB {

		tag "ALL"
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

	// Blast each protein chunk against the soft-masked genome
	// This is used to define targets for exhaustive exonerate alignments
	// has to run single-threaded due to bug in blast+ 2.5.0 (comes with Repeatmaster in Conda)
	process runBlastProteins {

		tag "Chunk ${chunk_name}"
		publishDir "${OUTDIR}/evidence/proteins/tblastn/chunks", mode: 'copy'

		input:
		file(protein_chunk) from fasta_prots
		file(blastdb_files) from blast_db_prots

		output:
		file(protein_blast_report) into ProteinBlastReport

		script:
		db_name = blastdb_files[0].baseName
		chunk_name = protein_chunk.getName().tokenize('.')[-2]
		protein_blast_report = "${protein_chunk.baseName}.blast"
		"""
			tblastn -db $db_name -query $protein_chunk -evalue $params.blast_evalue -outfmt "${params.blast_options}" > $protein_blast_report
		"""
	}

	// Parse Protein Blast output for exonerate processing
	process Blast2QueryTargetProts {

		tag "ALL"
	        publishDir "${OUTDIR}/evidence/proteins/tblastn/chunks", mode: 'copy'

		input:
		file(blast_reports) from ProteinBlastReport.collect()

		output:
		file(query2target_result_uniq_targets) into query2target_uniq_result_prots
	
		script:
		query_tag = Proteins.baseName
		query2target_result_uniq_targets = "${query_tag}.targets"
	
		"""
			cat $blast_reports > merged.txt
			blast2exonerate_targets.pl --infile merged.txt --max_intron_size $params.max_intron_size > $query2target_result_uniq_targets
		"""
	}

	// split Blast hits for parallel processing in exonerate
	query2target_uniq_result_prots
		.splitText(by: params.nexonerate, file: true)
		.combine(ProteinDB)
		.set{ query2target_chunk_prots }

	// Run Exonerate on the blast regions
	process runExonerateProts {

		tag "Chunk ${chunk_name}"
		// publishDir "${OUTDIR}/evidence/proteins/exonerate/chunks", mode: 'copy'

		input:
		set file(hits_chunk),file(protein_db),file(protein_db_index) from query2target_chunk_prots
		set file(genome),file(genome_faidx) from RMGenomeIndexProtein
	
		output:
		file(exonerate_chunk) into exonerate_result_prots
		file("merged.${chunk_name}.exonerate.out") into exonerate_raw_results
	
		script:
		query_tag = protein_db.baseName
		chunk_name = hits_chunk.getName().tokenize('.')[-2]
		exonerate_chunk = "${hits_chunk.baseName}.${query_tag}.exonerate.out"
		
		// get the protein fasta sequences, produce the exonerate command and genomic target interval fasta, run the whole thing,
		// merge it all down to one file and translate back to genomic coordinates
		// remove all the untracked intermediate files

		"""
			extractMatchTargetsFromIndex.pl --matches $hits_chunk --db $protein_db_index
			exonerate_from_blast_hits.pl --matches $hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --query_index $protein_db_index --analysis protein2genome --outfile commands.txt
			parallel -j ${task.cpus} < commands.txt
			cat *.exonerate.align | grep -v '#' | grep 'exonerate:protein2genome:local' > merged.${chunk_name}.exonerate.out
			exonerate_offset2genomic.pl --infile merged.${chunk_name}.exonerate.out --outfile $exonerate_chunk
			rm *.align
			rm *._target_.fa*
			rm *._query_.fa*
		"""
	}

	// merge the exonerate hits and create the hints
	process Exonerate2HintsProtein {

		tag "ALL"
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
			exonerate2gff.pl --infile all_chunks.out --source protein --outfile $exonerate_gff
		"""
	}

} // close protein loop

// --------------------------
// -------------------
// EST DATA PROCESSING
//-------------------
// --------------------------

if (params.ESTs != false ) {

	// create a cdbtools compatible database for ESTs
	process runIndexESTDB {

		tag "ALL"
		publishDir "${OUTDIR}/databases/cdbtools/ESTs", mode: 'copy'

		input:
		file (est_fa) from ests_index

		output:
		set file(est_fa),file(est_index) into EstDB

		script:
		est_index = est_fa.getName() + ".cidx"

		"""
			cdbfasta $est_fa
		"""
	}

	/*
	 * EST blasting
	*/

	// Blast each EST chunk against the nucleotide database
	process runBlastEst {

		tag "Chunk ${chunk_name}"
		// publishDir "${OUTDIR}/evidence/EST/blast/chunks"

		input:
		file(est_chunk) from fasta_ests
		file(blastdb_files) from blast_db_ests
		
		output:
		file(blast_report) into ESTBlastReport

		script:
		db_name = blastdb_files[0].baseName
		chunk_name = est_chunk.getName().tokenize('.')[-2]
		blast_report = "${est_chunk.baseName}.${db_name}.est.blast"

		"""
			blastn -db $db_name -evalue $params.blast_evalue -query $est_chunk -outfmt "${params.blast_options}" -num_threads ${task.cpus} > $blast_report
		"""
	}

	// Parse the EST Blast output
	process Blast2QueryTargetEST {

        	publishDir "${OUTDIR}/evidence/EST/blast", mode: 'copy'

		input:
		file(blast_report) from ESTBlastReport.collect()

		output:
		file(targets) into est_blast_targets

		script:
		query_tag = ESTs.baseName
		targets = "EST.blast.targets.txt"

		"""
			cat $blast_report >> merged.out
			blast2exonerate_targets.pl --infile merged.out --max_intron_size $params.max_intron_size > $targets
		"""

	}

	// Split EST targets and intersect with the Cdbtools index for fast target retrieval
	est_blast_targets
		.splitText( by: params.nexonerate , file: true )
		.combine(EstDB)
		.set { est_exonerate_chunk }

	// Run exonerate on the EST Blast chunks
	process runExonerateEST {

		tag "Chunk ${chunk_name}"
		// publishDir "${OUTDIR}/evidence/EST/exonerate/chunks", mode: 'copy'

		input:
		set file(est_hits_chunk),file(est_fa),file(est_db_index) from est_exonerate_chunk
		set file(genome),file(genome_index) from RMGenomeIndexEST	

		output:
		file(results) into exonerate_result_ests
	
		script:	
		chunk_name = est_hits_chunk.getName().tokenize('.')[-2]
		results = "EST.${chunk_name}.exonerate.out"	

		"""
			extractMatchTargetsFromIndex.pl --matches $est_hits_chunk --db $est_db_index
                        exonerate_from_blast_hits.pl --matches $est_hits_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --query_index $est_db_index --analysis est2genome --outfile commands.txt
                        parallel -j ${task.cpus} < commands.txt
                        cat *.exonerate.align |  grep "exonerate:est2genome" > merged_exonerate.out
                        exonerate_offset2genomic.pl --infile merged_exonerate.out --outfile $results
                        rm *.align
			rm *._target_.fa*
			rm *._query_.fa*
		"""
	}

	// Combine exonerate hits and generate hints
	process Exonerate2HintsEST {

		tag "ALL"
		publishDir "${OUTDIR}/evidence/EST/exonerate/", mode: 'copy'
	
		input:
		file(exonerate_result) from exonerate_result_ests.collect()
	
		output:
		file(exonerate_hints) into est_exonerate_hints
	
		script:
		exonerate_hints = "ESTs.exonerate.hints.gff"
			
		"""
			cat $exonerate_result >> merged.out
			exonerate2gff.pl --infile merged.out --source est --outfile $exonerate_hints
		"""
	}

} // close EST loop

// ++++++++++++++++++
// RNA-seq PROCESSING
// ++++++++++++++++++
if (params.reads != false ) {
	// ++++++++++++++++++
	// RNA-seq PROCESSING
	// ++++++++++++++++++

	// trim reads
	process runFastp {

		tag "${prefix}"
		publishDir "${OUTDIR}/evidence/rnaseq/fastp", mode: 'copy'

		input:
		set val(name), file(reads) from read_files_trimming

		output:
		file("*_trimmed.fastq.gz") into trimmed_reads
		set file(json),file(html) into trimmed_reads_qc

		script:
		prefix = reads[0].toString().split("_R1")[0]
		json = file(reads[0]).getBaseName() + ".fastp.json"
		html = file(reads[0]).getBaseName() + ".fastp.html"

		if (params.singleEnd) {
			left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
			"""
                       		fastp -i ${reads[0]} --out1 ${left} -w ${task.cpus} -j $json -h $html
	                """
		} else {
			left = file(reads[0]).getBaseName() + "_trimmed.fastq.gz"
			right = file(reads[1]).getBaseName() + "_trimmed.fastq.gz"
			"""
				fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html
			"""
		}
	}

	// Generate an alignment index from the genome sequence
	process runMakeHisatDB {

		tag "${prefix}"
		publishDir "${OUTDIR}/databases/HisatDB", mode: 'copy'

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
			hisat2-build $genome $dbName -p ${task.cpus}
			"""
		}	
	}

	/*
	 * STEP RNAseq.3 - Hisat2
	 */

	process runHisat2 {

		tag "${prefix}"
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
	process mergeHisatBams {

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
	process Hisat2Hints {
	
		// publishDir "${OUTDIR}/evidence/rnaseq/hints/chunks", mode: 'copy'

		input:
		file(bam) from Hisat2Hints
	
		output:
		file(hisat_hints) into rnaseq_hints
	
		script:
		hisat_hints = "RNAseq.hisat.hints.gff"

		"""
			bam2hints --intronsonly 0 -p 5 -s 'E' --in=$bam --out=$hisat_hints
		"""
	}

	// ----------------------------
	// run trinity de-novo assembly
	// ----------------------------

	if (params.trinity != false ) {

		process runTrinity {
	
			publishDir "${OUTDIR}/evidence/rnaseq/trinity", mode: 'copy'

			//scratch true 
	
			input:
			file(hisat_bam) from bam2trinity

			output:
			file "transcriptome_trinity/Trinity-GG.fasta" into trinity_transcripts, trinity_transcripts_2exonerate, trinity_to_index
	
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

		trinity_chunks = trinity_transcripts.splitFasta(by: params.nblast, file: true)

		process runBlastTrinity {
	
			tag "Chunk ${chunk_name}"
			// publishDir "${OUTDIR}/evidence/rnaseq/trinity/blast/chunks"
	
			input:
			file(query_fa) from trinity_chunks 
			file(blastdb) from blast_db_trinity
	
			output:
			file(blast_report) into TrinityBlastReport
	
			script: 

			db_name = blastdb[0].baseName
			chunk_name = query_fa.getName().tokenize('-')[-2]
			blast_report = "trinity.${chunk_name}.blast"

			"""
				blastn -db $db_name -query $query_fa -evalue $params.blast_evalue -outfmt "${params.blast_options}" -num_threads ${task.cpus} > $blast_report
			"""
		}

		/*
		 * STEP RNAseq.7 - Parse Blast Output
		 */

		process BlastTrinity2QueryTarget {
	
			publishDir "${OUTDIR}/evidence/rnaseq/trinity/exonerate", mode: 'copy'
	
			input:
			file(all_blast_results_trinity) from TrinityBlastReport.collectFile()
	
			output:
			file(trinity_targets) into query2target_trinity_uniq_result
	
			script:
			trinity_targets = "trinity.all.targets.txt"
			"""
	                        blast2exonerate_targets.pl --infile $all_blast_results_trinity --max_intron_size $params.max_intron_size > $trinity_targets
			"""
		} 	

		// generate an index for trinity transcripts with cdbtools
		process runTrinityIndex {

			tag "ALL"
			publishDir "${OUTDIR}/databases/trinity", mode: 'copy'

			input:
			file(trinity_fa) from trinity_to_index

			output:
			set file(trinity_fa),file(trinity_db_index) into TrinityDBIndex

			script:
			trinity_db_index = trinity_fa.getName() + ".cidx"
		
			"""
				cdbfasta $trinity_fa
			"""
	
		}

		// Split trinity targets and combine with trinity cdbtools index
		query2target_trinity_uniq_result
			.splitText(by: params.nexonerate, file: true)
			.combine(TrinityDBIndex)	
			.set{query2target_trinity_chunk}	

		/*
		 * STEP RNAseq.8 - Exonerate
		 */	
 
		process runExonerateTrinity {
	
			tag "Chunk ${chunk_name}"
			// publishDir "${OUTDIR}/evidence/rnaseq/trinity/exonerate/chunks", mode: 'copy'
	
			input:
			set file(hits_trinity_chunk),file(transcript_fa),file(transcript_db_index) from query2target_trinity_chunk
			set file(genome),file(genome_index) from RMGenomeIndexTrinity
		
			output:
			file(exonerate_out) into exonerate_result_trinity
	
			script:
			chunk_name = hits_trinity_chunk.getName().tokenize('.')[-2]
			exonerate_out = "RNAseq.Trinity.exonerate.${chunk_name}.out"
			
			"""
				extractMatchTargetsFromIndex.pl --matches $hits_trinity_chunk --db $transcript_db_index
	                        exonerate_from_blast_hits.pl --matches $hits_trinity_chunk --assembly_index $genome --max_intron_size $params.max_intron_size --query_index $transcript_db_index --analysis est2genome --outfile commands.txt
        	                parallel -j ${task.cpus} < commands.txt
                	        cat *.exonerate.align | grep -v '#' | grep 'exonerate:est2genome' > merged_exonerate.out
                        	exonerate_offset2genomic.pl --infile merged_exonerate.out --outfile $exonerate_out
	                        rm *.align
				rm *._target_.fa*
				rm *._query_.fa*	
			"""
		}

		/*
		 * STEP RNAseq.9 - Exonerate to Hints
		 */
 
		process Exonerate2HintsTrinity {	

			tag "ALL"
			publishDir "${OUTDIR}/evidence/rnaseq/trinity/exonerate/", mode: 'copy'

			input:
			file(exonerate_results) from  exonerate_result_trinity.collect()

			output:
			file(trinity_hints) into trinity_exonerate_hints

			script:
			trinity_hints = "RNAseq.trinity.hints.gff"
			"""
				cat $exonerate_results | grep -v '#' | grep 'exonerate:est2genome' > all_chunks.out
	                        exonerate2gff.pl --infile all_chunks.out --source trinity --outfile $trinity_hints
			"""
		}

	} // Close Trinity loop
	
} // Close RNAseq loop

/*
* RUN AUGUSTUS GENE PREDICTOR
*/

// get all available hints and merge into one file
process runMergeAllHints {

	tag "ALL"
	publishDir "${OUTDIR}/evidence/hints", mode: 'copy'

	input:
	file(protein_exonerate_hint) from prot_exonerate_hints.ifEmpty(false)
	file(rnaseq_hint) from rnaseq_hints.ifEmpty(false)
        file(est_exonerate_hint) from est_exonerate_hints.ifEmpty(false)
        file(trinity_exonerate_hint) from trinity_exonerate_hints.ifEmpty(false)

	output:
	file(merged_hints) into (mergedHints,mergedHintsSort)

	script:
	def file_list = ""
	if (protein_exonerate_hint != false || protein_exonerate_hint != "false" ) {
		file_list += " ${protein_exonerate_hint}"
	}
	if (rnaseq_hint  != false || rnaseq_hint != "false" ) {
		file_list += " ${rnaseq_hint}"
	}
	if (est_exonerate_hint != false || est_exonerate_hint != "false" ) {
		file_list += " ${est_exonerate_hint}"
	}
	if (trinity_exonerate_hint != false || trinity_exonerate_hint != "false" ) {
		file_list += " ${trinity_exonerate_hint}"
	}
	merged_hints = "merged.hints.gff"
	
	"""
		cat $file_list | grep -v "false" >> $merged_hints
	"""
}

process runHintsToBed {

	input:
	file(hints) from mergedHintsSort

	output:
	file(bed) into HintRegions

	script:

	bed = "regions.bed"

	"""
                grep -v "#" $hints | grep -v "false" | sort -k1,1 -k4,4n -k5,5n -t\$'\t' > hints.sorted
		gff2clusters.pl --infile hints.sorted --max_intron $params.max_intron_size > $bed
	"""
}

// execute Augustus
/*
 * STEP Augustus.1 - Genome Annotation
 */
// Run against each repeatmasked chunk of the assembly
process runAugustus {

	tag "Chunk ${chunk_name}"
	//publishDir "${OUTDIR}/annotation/augustus/chunks"

        when:
        params.augustus != false

	input:
	file(regions) from HintRegions.collect()
	file(hints) from mergedHints.collect()
	file(genome_chunk) from GenomeChunksAugustus

	output:
	file(augustus_result) into augustus_out_gff
	
	script:
	chunk_name = genome_chunk.getName().tokenize(".")[-2]
	augustus_result = "augustus.${chunk_name}.out.gff"
	genome_fai = genome_chunk.getName() + ".fai"

	"""
		samtools faidx $genome_chunk
		fastaexplode -f $genome_chunk -d . 
		augustus_from_regions.pl --genome_fai $genome_fai --model $params.model --utr $params.UTR --isof $params.isof --aug_conf $AUG_CONF --hints $hints --bed $regions > commands.txt	
		parallel -j ${task.cpus} < commands.txt
		touch dummy.augustus.gff
		cat *augustus.gff > $augustus_result
	"""
}

process runMergeAugustusGff {

	tag "ALL"
	publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'
	
	input:
	file(augustus_gffs) from augustus_out_gff.collect()

	output:
	file(augustus_merged_gff) into (augustus_2gff3, augustus_2prots)

	script:
	augustus_merged_gff = "augustus.merged.out.gff"
	
	"""	
		cat $augustus_gffs >> merged.gff
		create_gff_ids.pl --gff merged.gff > $augustus_merged_gff
	"""
}

process runAugustus2Protein {

	tag "ALL"
        publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'

	input:
	file augustus2parse from augustus_2prots
	
	output:
	file(augustus_prot_fa) into (augustus_proteins, augustus_prots2annie, augustus_prots2interpro)
	
	script:
	augustus_prot_fa = "augustus.proteins.fa"

	"""
		getAnnoFasta.pl $augustus2parse
		cat *.aa > $augustus_prot_fa	
	"""
}

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
