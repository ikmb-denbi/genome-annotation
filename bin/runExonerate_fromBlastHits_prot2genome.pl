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
my $genome_file = $ARGV[2];		# GENOME FILE
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

# CREATE DIRECTORY WITH ONE FILE PER SCAFFOLD
if (-e $genome_basename and -d $genome_basename) {

} else {
	system("mkdir $genome_basename");

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

	my $query_clean = $query;
	$query_clean =~ s/\|/_/g ; # Strip characters that will break bash
	my $fa_file = "$query_clean.fa" ;
	
	#printf "cdbyank $queries_fasta -a '$query' > $fa_file \n"; 
	
	system("cdbyank $queries_fasta -a '$query' > $fa_file");

	printf "exonerate --model protein2genome --softmaskquery yes --softmasktarget yes --bestn 1 --minintron 20 --maxintron 20000  --showalignment false --showtargetgff true $fa_file $genome_basename/$target.fasta >> $query_clean.exonerate.out\n";
}

close MATCHES;
exit;
