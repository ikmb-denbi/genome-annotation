#!/usr/bin/env nextflow
/*
========================================================================================
                    IKMB - de.NBI | Genome Annotation Pipeline
========================================================================================
 Genome Annotation Pipeline. Started 2018-10-17.
 #### Homepage / Documentation
 https://git.ikmb.uni-kiel.de/m.torres/NF-hints.git
 #### Authors
 MTorres m.torres <m.torres@ikmb.uni-kiel.de> - https://git.ikmb.uni-kiel.de/m.torres>
 MHoeppner m.hoeppner <m.hoeppner@ikmb.uni-kiel.de> 
----------------------------------------------------------------------------------------
*/

// Make sure the Nextflow version is current enough
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nextflow_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please use a more recent version of Nextflow!\n" +
              "============================================================"
}

def helpMessage() {
  log.info"""
  =================================================================
   IKMB - de.NBI | Genome Annotation Pipeline | v${nextflow.manifest.version}
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
    --gth		Run GenomeThreader to produce hints from protein file [ true (default) | false ]
    --augustus		Run Augustus to predict genes [ true (default) | false ]
    --funAnnot		Run functional annotation using Annie [ true (default) | false ]
 	
    Programs parameters:
    --species		Species database for RepeatMasker [ default = 'mammal' ]
    --rm_lib		Additional repeatmasker library in FASTA format [ default = 'false' ]
    --model		Species model for Augustus [ default = 'human' ]
    --UTR		Allow Augustus to predict UTRs (results are not optimal and takes much longer) [ 'on' | 'off' (default) ]
    --iso		Allow Augustus to predict multiple isoforms  (results are not optimal and takes much longer) [ 'true' | 'false' (default) ]
    --augCfg		Location of augustus configuration file [ default = 'bin/augustus_default.cfg' ]
    --uniprot		Fasta file with Uniprot proteins for functional annotation [ default = '/bin/Eumetazoa_UniProt_reviewed_evidence.fa' ]
 	
    How to split programs:
    --nblast		Chunks (# of sequences) to divide Blast jobs [ default = 500 ]
    --nexonerate	Chunks (# of blast hits) to divide Exonerate jobs [ default = 200 ]
    --nrepeats		Chunks (# of scaffolds) to divide RepeatMasker and Augustus jobs [ default = 30 ]
    --ninterpro		Chunks (# of sequences) to divide InterPro jobs [ default = 200 ]

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
	if(params.gth == false) {
		println "Will run only Exonerate on Protein file"
	} else {
                println "Will run Exonerate and GenomeThreader on Protein file"
	}
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
}

// Make it fail if basic requirements are unmet
if (!binding.variables.containsKey("Proteins") && !binding.variables.containsKey("ESTs") && params.reads == false) {
	exit 1, "At least one type of input data must be specified (--proteins, --ESTs, --reads)"
}

if (params.trinity == true && params.reads == false) {
	exit 1, "Cannot run Trinity de-novo assembly without RNA-seq reads (specify both --reads and --trinity)"
}

// give this run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

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
// --- align proteins with genomethreader
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
        Channel
        .fromPath(Proteins)
        .splitFasta(by: params.nblast, file: true)
        .set { fasta_prots }
	Channel
	.fromPath(Proteins)
	.set { index_prots }
} else {
        fasta_prots = Channel.from(false)
        prot_exonerate_hints = Channel.create()
	index_prot = Channel.from(false)
}

// if ESTs are provided
if (params.ESTs) {
        Channel
                .fromPath(ESTs)
                .splitFasta(by: params.nblast, file: true)
                .set {fasta_ests}
	Channel
		.fromPath(ESTs)
		.set { ests_index }
} else {
        fasta_ests = Channel.from(false)
	ests_index = Channel.from(false)
        est_exonerate_hints = Channel.create()
}

// if GenomeThreader should be run
if (params.gth == false ) {
        fasta_prots_gth = Channel.from(false)
        gth_protein_hints = Channel.create()
} else {
        Channel
        .fromPath(Proteins)
        .splitFasta(by: params.nblast, file: true)
        .set {fasta_prots_gth}
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
		.into {read_files_trimming }
} else {
	// this channel should emmit two items...
	read_files_trimming = Channel.from([false,false])
	inputMakeHisatdb = Channel.from(false)
	rnaseq_hints = Channel.from(false)
	trinity_exonerate_hints = Channel.from(false)
}

// Header log info
log.info "========================================="
log.info "IKMB Genome Annotation Pipeline v${workflow.manifest.version}"
log.info "Genome assembly: 		${params.genome}"
log.info "-----------------------------------------"
log.info "Evidences:"
log.info "Proteins:			${params.proteins}"
log.info "ESTs:				${params.ESTs}"
log.info "RNA-seq:			${params.reads}"
log.info "-----------------------------------------"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="

// ---------------------------
// RUN REPEATMASKER
//----------------------------

// generate a soft-masked sequence for each assembly chunk
process runRepeatMasker {

	tag "Chunk ${chunk_name}"
	publishDir "${OUTDIR}/repeatmasker/chunks"

	input: 
	file(genome_fa) from FastaRM

	output:
	file(genome_rm) into RMFastaChunks

	script:
	chunk_name = genome_fa.getName().tokenize('.')[-2]
	// provide a custom repeat mask database
	options = ""
	if (params.rm_lib) {
		options = "-lib ${RM_LIB}"
	}
	genome_rm = genome_fa + ".masked"
	
	"""
		RepeatMasker -species ${params.species} -gff -xsmall  $options -q -pa ${task.cpus} $genome_fa	
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

	script:
	
	masked_genome = "${Genome.baseName}.rm.fa"
	"""
		cat $genome_chunks >> merged.fa
		fastasort -f merged.fa > $masked_genome
		rm merged.fa
	"""	
}
GenomeChunksAugustus = RMtoPartition
	.splitFasta(by: params.nrepeats, file: true)

// Split repeat-masked genome into individual sequences for efficient exonerate alignments
process runExplodeGenome {

	tag "ALL"
	publishDir "${OUTDIR}/databases/cdbtools/genome", mode: 'copy'

	input:
	file(genome_fa) from RMtoSplit

	output:
	file(genome_chunk_dir) into (GenomeChunksProtein,GenomeChunksEst,GenomeChunksTrinity)

	script:
	genome_chunk_dir = "genome"

	"""
		mkdir -p $genome_chunk_dir
		fastaexplode -f $genome_fa -d $genome_chunk_dir
	"""
}

// Turn genome into a masked blast database
// Generates a dust mask from softmasked genome sequence
process runMakeBlastDB {
	
	tag "ALL"
	publishDir "${OUTDIR}/databases/blast/", mode: 'copy'

	input:
	file(genome_fa) from RMtoBlastDB

	output:
	file("${dbName}.n*") into (blast_db_prots, blast_db_ests, blast_db_trinity)
	file(db_mask) into BlastDBMask

	script:
	dbName = genome_fa.baseName
	db_mask = "${dbName}.asnb"
	
	"""
		convert2blastmask -in $genome_fa -parse_seqids -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out $db_mask
		makeblastdb -in $genome_fa -parse_seqids -mask_data $db_mask -dbtype nucl -out $dbName 
	"""
}

if (params.proteins != false) {
	// ---------------------
	// PROTEIN DATA PROCESSING
	// ---------------------

	// ----------------------------
	// Protein BLAST against genome
	// ----------------------------

	// create a cdbtools compatible  index
	// we need this to do very focused exonerate searches later
	process runIndexProteinDB {

		tag "ALL"
		publishDir "${OUTDIR}/databases/cdbtools/proteins"

		input:
		file(fasta) from index_prots

		output:
		set file(fasta),file(protein_index) into ProteinDB

		when: 
		params.proteins != false

		script:
		protein_index = fasta.getName()+ ".cidx"

		"""
			cdbfasta $fasta 
		"""
	}

	// Blast each protein chunk against the soft-masked genome
	// This is used to define targets for exhaustive exonerate alignments
	process runBlastProteins {

		tag "Chunk ${chunk_name}"
		publishDir "${OUTDIR}/evidence/proteins/tblastn/chunks"

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
			tblastn -db $db_name -query $protein_chunk -db_soft_mask 40 -max_target_seqs 1 -outfmt 6 > $protein_blast_report
		"""
	}

	// Parse Protein Blast output for exonerate processing
	process Blast2QueryTargetProts {

		tag "ALL"
	        publishDir "${OUTDIR}/evidence/proteins/tblastn/chunk/"

		input:
		file(blast_reports) from ProteinBlastReport.collect()

		output:
		file(query2target_result_uniq_targets) into query2target_uniq_result_prots
	
		script:
		query_tag = Proteins.baseName
		query2target_result_uniq_targets = "${query_tag}.targets"
	
		"""
			cat $blast_reports > merged.txt
			BlastOutput2QueryTarget.pl merged.txt 1e-5 query2target_result
			sort query2target_result | uniq > $query2target_result_uniq_targets
		"""
	}

	// split Blast hits for parallel processing in exonerate
	query2target_uniq_result_prots
		.splitText(by: params.nexonerate, file: true)
		.combine(ProteinDB)
		.set{ query2target_chunk_prots }

	// Run Exonerate on the blast regions
	process runExonerateProts {

		tag "${query_tag}"
		publishDir "${OUTDIR}/evidence/proteins/exonerate/chunks"
	
		input:
		set file(hits_chunk),file(protein_db),file(protein_db_index) from query2target_chunk_prots
		file(genome_root) from GenomeChunksProtein
	
		output:
		file(exonerate_chunk) into exonerate_result_prots
	
		script:
		query_tag = protein_db.baseName
		exonerate_chunk = "${hits_chunk.baseName}.${query_tag}.exonerate.out"
	
		"""
			runExonerate_fromBlastHits_prot2genome.pl $hits_chunk $protein_db_index $genome_root >> commands.txt
			parallel -j ${task.cpus} < commands.txt
			cat *.exonerate.out | grep -v '#' | grep 'exonerate:protein2genome:local' > $exonerate_chunk
		"""
	}

	// merge the exonerate hits and create the hints
	process runMergeExonerateHints {

		publishDir "${OUTDIR}/evidence/proteins/exonerate/", mode: 'copy'

		input:
		file(chunks) from exonerate_result_prots.collect()

		output:
		file(exonerate_hints) into prot_exonerate_hints

		script:
		query_tag = Proteins.baseName
		exonerate_gff = "proteins.exonerate.${query_tag}.hints.gff"
		"""
			cat $chunks > all_chunks.out
			Exonerate2GFF_protein.pl all_chunks.out $exonerate_gff
		"""
	}

	// ------------------------------------
	// GenomeThreader hints generation
	// ------------------------------------

	// Run genome threader for protein chunks if requested
	process runGenomeThreaderProteins {

		tag "Chunk ${chunk_name}"
		publishDir "${OUTDIR}/evidence/proteins/gth/chunks/"

		input:
		file(protein_chunk) from fasta_prots_gth

		output:
		file(gth_chunk) into ProteinGTHChunk
	
		when:
		params.gth 

		script:
		chunk_name = protein_chunk.getName().tokenize('.')[-2]
		gth_chunk = "${protein_chunk.getName()}.gth"

		"""
			gth -genomic $Genome -protein $protein_chunk -gff3out -intermediate -o $gth_chunk
		"""	
	}

	// convert gth hits into hints
	process GenomeThreader2HintsProts {

		tag "${chunk_name}"
	        publishDir "${OUTDIR}/evidence/proteins/gth/chunks"

		input:
		file(not_clean_gth) from ProteinGTHChunk
	
		output:
		file(gth_hints) into ProteinGTHChunkHint
	
		script:
		gth_hints = not_clean_gth.baseName + ".clean.gth"
		chunk_name = not_clean_gth.getName().tokenize('.')[-2]

		"""
			gt gff3 -addintrons yes -setsource gth -tidy yes -addids no $not_clean_gth > not_clean_gth_wIntrons
			grep -v '#' not_clean_gth_wIntrons > no_hash_gth
			GTH_rename_splitoutput.pl no_hash_gth > clean_gth
			grep -e 'CDS' -e 'exon' -e 'intron' clean_gth | perl -ple 's/Parent=/grp=/' | perl -ple 's/(.*)\$/\$1;src=P;pri=3/' | perl -ple 's/CDS/CDSpart/' | perl -ple 's/intron/intronpart/' | perl -ple 's/exon/exonpart/' > $gth_hints
		
		"""
	}

	// merge gth hints
	process GenomeThreaderMergeHints {

		tag "ALL"
        	publishDir "${OUTDIR}/evidence/proteins/gth/", mode: 'copy'

		input:
		file(gth_hint_chunks) from ProteinGTHChunkHint.collect()

		output:
		file(merged_gth_hints) into gth_protein_hints

		script:
        	query_tag = Proteins.baseName
		merged_gth_hints = "proteins.gth.${query_tag}.hints.gff"

		"""
			cat $gth_hint_chunks >> $merged_gth_hints
		"""
	}

} // close protein loop

if (params.ESTs != false) {

	// --------------------------
	// -------------------
	// EST DATA PROCESSING
	//-------------------
	// --------------------------

	// create a cdbtools compatible database for ESTs
	process runIndexESTDB {

		tag "ALL"
		publishDir "${OUTDIR}/databases/cdbtools/ESTs", mode: 'copy'

		input:
		file (est_fa) from ests_index

		output:
		set file(est_fa),file(est_index) into EstDB

		when: 
		params.ESTs != false

		script:
		est_index = est_fa.getName() + ".cdix"

		"""
			cdbfasta $est_fa
		"""
	}

	/*
	 * EST blasting
	*/

	// Blast each EST chunk against the nucleotide database
	process runBlastEst {

		tag "Chunk: #{chunk_name}"
		publishDir "${OUTDIR}/evidence/EST/blast/chunks"

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
			blastn -db $db_name -db_soft_mask 40 -query $query_fa_ests -max_target_seqs 1 -outfmt 6 -num_threads ${task.cpus} > $blast_report
		"""
	}

	// Parse the EST Blast output
	process Blast2QueryTargetEST {

        	publishDir "${OUTDIR}/evidence/EST/blast", mode: 'copy'

		input:
		file(blast_report) from ESTBlastReport.collectFile()

		output:
		file(query2target_result_uniq_targets) into query2target_uniq_result_ests

		script:
		query_tag = ESTs.baseName
		query2target_result_uniq_targets = "EST.blast.targets.txt"

		"""
			BlastOutput2QueryTarget.pl $blast_report 1e-5 query2target_result
			sort query2target_result | uniq > $query2target_result_uniq_targets
		"""

	}

	// Split EST targets and intersect with the Cdbtools index for fast target retrieval
	query2target_uniq_result_ests
		.splitText(by: params.nexonerate, file: true)
		.combine(EstDB)
		.set{query2target_chunk_ests}

	// Run exonerate on the EST Blast chunks
	process runExonerateEST {

		tag "${chunk_name}"
		publishDir "${OUTDIR}/evidence/EST/exonerate/chunks"
	
		input:
		set file(hits_chunk),file(est_fa),file(est_db_index) from query2target_chunk_ests
		file(genome_base_dir) from GenomeChunksEst	

		output:
		file(results) into exonerate_result_ests
	
		script:	
		query_tag = ESTs.baseName
		chunk_name = hits_chunk.getName().tokenize('.')[-1]
		results = "${hits_chunk.getName()}.${query_tag}.exonerate.out"
	
		"""
			runExonerate_fromBlastHits_est2genome.pl $hits_chunk $est_db_index $genome_base_dir >> commands.txt
			parallel -j ${task.cpus} < commands.txt
			cat *.exonerate.out >> $results
		"""
	}

	// Combine exonerate hits and generate hints
	process Exonerate2HintsEST {

		tag "${query_tag}"
		publishDir "${OUTDIR}/evidence/EST/exonerate/", mode: 'copy'
	
		input:
		file(exonerate_result) from exonerate_result_ests.collectFile()
	
		output:
		file(exonerate_gff) into est_exonerate_hints
	
		script:
		query_tag = ESTs.baseName
		exonerate_hints = "ESTs.exonerate.${query_tag}.hints.gff"
			
		"""
			grep -v '#' $exonerate_result_ests | grep 'exonerate:est2genome' > exonerate_gff_lines
			Exonerate2GFF_EST.pl exonerate_gff_lines $exonerate_hints
		"""
	}

} // close EST loop

if (params.reads != false) {
	// ++++++++++++++++++
	// RNA-seq PROCESSING
	// ++++++++++++++++++

	/*
	 * Create a channel for input read files
	 */
 
	/*
	 * STEP RNAseq.1 - Trimgalore
	 */

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
		file "*accepted_hits.bam" into accepted_hits2merge , bam2trinity
	
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

	process mergeHisatBams {

		publishDir "${OUTDIR}/evidence/ranseq/Hisat2", mode: 'copy'

		input:
		file hisat_bams from accepted_hits2merge.collect()

		output:
		file(bam) into  Hisat2Hints

		script:
		bam = "hisat2.merged.bam"
		avail_ram_per_core = (task.memory/${task.cpus}).toGiga()-1
	
		"""
			samtools merge - $hisat_bams | samtools sort -@ ${task.cputs} -m${avail_ram_per_core}G - > $bam
		"""
	}

	/*
	 * STEP RNAseq.4 - Hisat2 into Hints
	 */	
	process Hisat2Hints {
	
		tag "${prefix}"
		publishDir "${OUTDIR}/evidence/rnaseq/hints/chunks", mode: 'copy'

		input:
		file(bam) from Hisat2Hints
	
		output:
		file(hisat_hints) into rnaseq_hints
	
		script:
		hisat_hints = "rnaseq.hisat.hints.gff"

		"""
			bam2hints --intronsonly 0 -p 5 -s 'E' --in=$accepted_hits2hints --out=${hisat_hints_}
		"""
	}

	// ----------------------------
	// run trinity de-novo assembly
	// ----------------------------

	process runTrinity {
	
		publishDir "${OUTDIR}/evidence/rnaseq/trinity", mode: 'copy'

		scratch true 
	
		input:
		file(hisat_bam) from bam2trinity.collect()

		output:
		file "transcriptome_trinity/Trinity-GG.fasta" into trinity_transcripts, trinity_transcripts_2exonerate, trinity_to_index
	
		script:

		trinity_option = ( params.rnaseq_stranded == true ) ? "--SS_lib_type RF" : ""

		"""
			Trinity --genome_guided_bam $hisat_bam \
			--genome_guided_max_intron 10000 \
			--CPU ${task.cpus} \
			--max_memory ${task.memory.toGiga()-1}G \
			--output transcriptome_trinity \
			$trinity_option
		"""
	}

	trinity_chunks = trinity_transcripts.splitFasta(by: params.nblast, file: true)

	process runBlastTrinity {
	
		tag "Chunk ${chunk_name}"
		publishDir "${OUTDIR}/evidence/rnaseq/trinity/blast/chunks"

		input:
		file(query_fa) from trinity_chunks 
		file(blastdb) from blast_db_trinity
	
		output:
		file(blast_report) into TrinityBlastReport
	
		script: 

		db_name = blastdb_nhr.baseName
		chunk_name = query_fa.getName().tokenize('-')[-2]
		blast_report = "trinity.${chunk_name}.blast"

		"""
			blastn -db $db_name -query $query_fa -db_soft_mask 40 -max_target_seqs 1 -outfmt 6 > blast_report
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
			BlastOutput2QueryTarget.pl $all_blast_results_trinity 1e-5 query2target_trinity_result
			sort query2target_trinity_result | uniq > $trinity_targets
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
		trinity_db_index = trinity_fa.getName() + ".cdix"
	
		"""
			cdbtools $trinity_fa
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
		publishDir "${OUTDIR}/evidence/rnaseq/trinity/exonerate/chunks", mode: 'copy'
	
		input:
		set file(hits_trinity_chunk),file(transcript_fa),file(transcript_db_index) from query2target_trinity_chunk
		file(genome_base_dir) from GenomeChunksTrinity
	
		output:
		file(exonerate_out) into exonerate_result_trinity
	
		script:
		chunk_name = hits_trinity_chunk.getName().tokenize('.')[-2]
		exonerate_out = "RNAseq.Trinity.exonerate.${chunk_name}.out"
		
		"""
			runExonerate_fromBlastHits_est2genome.pl $hits_trinity_chunk $transcript_db_index $genome_base_dir >> commands.txt
			parallel -j ${task.cpus} < cat commands.txt
			cat *.exonerate.out >> $exonerate_out
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
			cat $exonerate_results | grep -v '#' | grep 'exonerate:est2genome' > exonerate_gff_lines
			Exonerate2GFF_trinity.pl exonerate_gff_lines $trinity_hints
		"""
	}

} // Close RNAseq loop

/*
* RUN AUGUSTUS GENE PREDICTOR
*/

// get all available hints and merge into one file
process runMergeAllHints {

	tag "ALL"
	publishDir "${OUTDIR}/evidence/hints", mode: 'copy'

	input:
	file(protein_exonerate_hint) from prot_exonerate_hints.ifEmpty()
	file(protein_gth_hint) from gth_protein_hints.ifEmpty()
	file(est_exonerate_hint) from est_exonerate_hints.ifEmpty()
	file(trinity_exonerate_hint) from trinity_exonerate_hints.ifEmpty()
	file(rnaseq_hint) from rnaseq_hints.ifEmpty()

	output:
	
	file(merged_hints) into mergedHints

	script:
	merged_hints "merged_all_hints.gff"

	"""
		cat $rnaseq_hint $protein_exonerate_hint $protein_gth_hint $est_exonerate_hint $trinity_exonerate_hint >> $merged_hints
	"""
}

// execute Augustus
/*
 * STEP Augustus.1 - Genome Annotation
 */
// Run against each repeatmasked chunk of the assembly
process runAugustus {

	tag "Chunk ${chunk_name}"
	publishDir "${OUTDIR}/annotation/augustus/chunks"

        when:
        params.augustus != false

	input:

	file(hints) from mergedHints
	file(genome_chunk) from GenomeChunksAugustus

	output:
	fil(augustus_out) into augustus_out
	
	script:
	chunk_name = genome_chunk.getName().tokenize(".")[-2]
	augustus_out = "augustus.${chunk_name}.out.gff"
	"""
		augustus --species=$params.model --UTR=$params.UTR --alternatives-from-evidence=$params.isof --extrinsicCfgFile=$AUG_CONF --hintsfile=$hints $Genome > $augustus_out
	"""
}

process runMergeAugustusGff {

	tag "ALL"
	publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'
	
	input:
	file(augustus_gffs) from augustus_out_gff.collect()

	output:
	file(augustus_merged_gff) into augustus_2gff3, augustus_2prots

	script:
	augustus_merged_gff = "augustus.merged.out.gff"
	
	"""
		cat $augustus_gffs > $augustus_merged_gff
	"""
}

process runAugustus2Protein {

	tag "ALL"
        publishDir "${OUTDIR}/annotation/augustus", mode: 'copy'

	when:
        params.augustus != false

	input:
	file augustus2parse from augustus_2prots
	
	output:
	file(augustus_prot_fa) into augustus_proteins, augustus_prots2annie, augustus_prots2interpro
	
	script:
	augustus_prot_fa = "augustus.proteins.fa"

	"""
		getAnnoFasta.pl $augustus2parse
		cat *.aa > $augustus_prot_fa	
	"""
}

