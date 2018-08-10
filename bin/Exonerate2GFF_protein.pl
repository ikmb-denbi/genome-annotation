#!/usr/bin/env perl

# By M. Torres on 26.10.2017

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 2)

{
    print "Usage of Exonerate2GFF.pl\n\n";
    print "perl Exonerate2GFF_both.pl <exonerate_output> <output_file>\n";
    print "where <exonerate_output> is the input file with exonerate results,\n";
    print " <output_file> is the output file for this script.\n";
    print "For example, >perl Exonerate2GFF.pl exonerate_output.out exonerate.gff\n";
    exit;
}

my $exonerate_out			= $ARGV[0];		# EXONERATE OUTPUT FILE
my $output_file				= $ARGV[1]; 	# OUTPUT FILE
my $line;
my @temp;
my $GeneID;									# ID OF THAT GENE
# COLUMNS FROM EXONERATE GFF:
my $Chrom;
my $method;
my $feature;
my $start;
my $end;
my $score;
my $strand;
my $frame;
my $comment;


# OPEN FILE WITH BLAST RESULTS:  
open(EXONERATE,"$exonerate_out") || die "ERROR: cannot open $exonerate_out\n";

# OPEN FILE TO WRITE OUTPUT ON:
open(OUTFILE,">$output_file") || die "ERROR: cannot open $output_file\n";

while (<EXONERATE>) {
	$line = $_;
	chomp $line;
	
	@temp =split(/\t+/,$line);
	$Chrom = $temp[0];
	$method = $temp[1];
	$feature = $temp[2];
	$start = $temp[3];
	$end = $temp[4];
	$score = $temp[5];
	$strand = $temp[6];
	$frame = $temp[7];
	$comment = $temp[8];
	
	($method) = ($method =~/exonerate:(\w+)/);
	if ($feature eq "gene") {
		($GeneID) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)\s;\s/);
	} elsif ($feature eq "exon") {
		print OUTFILE $Chrom."\t".$method."\texonpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=P;pri=3\n";
	} elsif ($feature eq "cds") {
		print OUTFILE $Chrom."\t".$method."\tCDSpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=P;pri=3\n";
	} elsif ($feature eq "intron") {
		print OUTFILE $Chrom."\t".$method."\tintronpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=P;pri=3\n";
	} else {
		next;
	}
}
close (EXONERATE);
close (OUTFILE);

exit