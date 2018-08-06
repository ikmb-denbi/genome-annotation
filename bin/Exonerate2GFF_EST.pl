#!/usr/bin/env perl

# By M. Torres on 07.11.2017

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 3)

{
    print "Usage of Exonerate2GFF_EST.pl\n\n";
    print "perl Exonerate2GFF_EST.pl <exonerate_output> <'var'|'no_var'> <output_file>\n";
    print "where <exonerate_output> is the input file with exonerate results,\n";
    print " <'var'|'no_var'> if there are variants or not in the queries,\n";
    print " <output_file> is the output file for this script.\n";
    print "For example, >perl Exonerate2GFF_EST.pl exonerate_output.out exonerate.gff\n";
    exit;
}

my $exonerate_out			= $ARGV[0];		# EXONERATE OUTPUT FILE
my $is_variant				= $ARGV[1];		# VARIANTS ['var' | 'no_var']
my $output_file				= $ARGV[2]; 	# OUTPUT FILE
my $line;
my @temp;
my $GeneID;									# ID OF THAT GENE
my $variant;								# VARIANT [-RA | -RB | ...]
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
	
	if ($is_variant eq "var") {
		if ($feature eq "gene") {
			($GeneID, $variant) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)-(\w+)\s;\s/);
			print OUTFILE $Chrom."\t".$method."\texpressed_sequence_match\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tID=".$GeneID.";Name=".$GeneID."\n";
		} elsif ($feature eq "exon") {
			print OUTFILE $Chrom."\t".$method."\tmatch_part\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tParent=".$GeneID."\n";
		} else {
			next;
		}
	} elsif ($is_variant eq "no_var") {
		if ($feature eq "gene") {
			($GeneID) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)\s;\s/);
			print OUTFILE $Chrom."\t".$method."\texpressed_sequence_match\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tID=".$GeneID.";Name=".$GeneID."\n";
		} elsif ($feature eq "exon") {
			print OUTFILE $Chrom."\t".$method."\tmatch_part\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tParent=".$GeneID."\n";
		} else {
			next;
		}
	}
}
close (EXONERATE);
close (OUTFILE);

exit