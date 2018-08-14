#!/bin/env perl

# By M. Torres on 24.10.2017

use strict;
use warnings;
use Bio::SeqIO;

my $num_args               = $#ARGV + 1;
if ($num_args != 3)
{
    print "Usage of runExonerate_fromBlastHits.pl\n\n";
    print "perl runExonerate_fromBlastHits.pl <parsed_blast_output> <queries_file> <genome_file>\n";
    print "where <parsed_blast_output> is the input file with unique query-scaffold combinations from blast results,\n";
    print "	<queries_file> is the fasta file with all the possible query sequences\n";
    print " <genome_file> is the complete genome file\n";
    print "For example, >perl runExonerate_fromBlastHits.pl query_target_uniq.txt queries.fa genome.fa\n";
    exit;
}

my $blast_output_parsed = $ARGV[0];	# FILE WITH UNIQUE QUERY - SCAFFOLD COMBINATIONS
my $queries_fasta = $ARGV[1];		# FASTA FILE WITH ALL QUERIES
my $genome_file = $ARGV[2];			# GENOME FILE
my $genome;							# GENOME FILE BASENAME
my $line;
my @temp;
my $query;							# PROTEIN QUERY OF THAT HIT
my $target;							# TARGET SCAFFOLD OF THAT HIT
my $scaffold;						# THAT SCAFFOLD NAME
my $seqs;							# THAT SCAFFOLD SEQUENCE
my $genome_basename;
my $seqio_obj;						# SEQIO OBJECT WITH ALL QUERIES
my $seq_obj;						# SEQIO SEQUENCE OBJECT WITH THAT QUERY
my $seqout;						# SEQIO OUTFILE OBJECT


open(GENOME, "<$genome_file") || die "ERROR: cannot open $genome_file";

$genome_basename = $genome_file;
$genome_basename =~ s/\..*$//;
print STDOUT "$genome_basename\n";

# CREATE DIRECTORY WITH ONE FILE PER SCAFFOLD
if (-e $genome_basename and -d $genome_basename) {
	print STDOUT "Directory $genome_basename already exists\n";
} else {
	system("mkdir $genome_basename");
	print STDOUT "Created directory $genome_basename\n";

	while (<GENOME>) {
		$line = $_;
		chomp $line;
		
    	# WRITE THE HEADER:
   		if ($line=~/^>/){
			($scaffold) = ($line=~ /^>(\S+)/);
			open (SCFILE, ">$genome_basename/$scaffold.fasta");
			print STDOUT "Created file $scaffold.fasta in directory $genome_basename\n";
			print SCFILE ">$scaffold\n";
			next;
			
		# WRITE THE SEQUENCE:
		} else {
			$seqs = $line;
			print SCFILE "$seqs\n";
		}
		next;
	}
	close GENOME;
}

# RUN EXONERATE FOR EACH PROTEIN-SCAFFOLD MATCH:
open(MATCHES, "<$blast_output_parsed") || die "ERROR: cannot open $blast_output_parsed";


# READ EACH HIT ENTRY: 
while (<MATCHES>){
	$line = $_;
	chomp $line;
	@temp = split(/\t+/,$line);
	$query = $temp[0];
    $target = $temp[1];
	
	# MAKE A FILE FOR THE SEQUENCE OF THIS QUERY:
	$seqio_obj = Bio::SeqIO->new('-file' => $queries_fasta, '-format' => 'fasta' );	
	while ($seq_obj = $seqio_obj->next_seq ) {
    	if ($seq_obj->display_id() eq $query) {
    	# CREATE THE FASTA FILE FOR THIS QUERY:
    	$seqout = Bio::SeqIO->new(-file => ">$query.fa", -format => "fasta");
    	# PRINT THE QUERY SEQUENCE:
    	print $seqout->write_seq($seq_obj);
    	}
	}
	print "Exonerating $query against $target ...\n";
	#system("echo \"$query\" | cdbyank $cdb_prot_index > $query.fa");
	system("exonerate --model est2genome --softmaskquery yes --softmasktarget yes --bestn 1 --minintron 20 --maxintron 20000  --showalignment false --showtargetgff true $query.fa $genome_basename/$target.fasta >> $query_exonerate.out");
	system("rm $query.fa");
}

close MATCHES;
exit;
