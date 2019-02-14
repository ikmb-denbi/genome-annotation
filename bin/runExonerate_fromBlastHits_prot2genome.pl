#!/bin/env perl

# By M. Torres on 24.10.2017

use strict;
use warnings;

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
my $queries_fasta = $ARGV[1];		# FASTA FILE CDBTOOLS INDEX
my $genome_basename = $ARGV[2];		# PATH TO GENOME SEQUENCES
my $line;
my @temp;
my $query;							# PROTEIN QUERY OF THAT HIT
my $target;							# TARGET SCAFFOLD OF THAT HIT

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

	unless (-e $fa_file) {	
		system("cdbyank $queries_fasta -a '$query' > $fa_file");
	}

	printf "exonerate --model protein2genome --softmasktarget yes --bestn 1 --minintron 20 --maxintron 50000  --showalignment false --showtargetgff true $fa_file $genome_basename/$target.fa >> $query_clean.exonerate.out\n";
}

close MATCHES;
exit;
