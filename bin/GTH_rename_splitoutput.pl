#!/usr/bin/env perl

# By M. Torres on 14.06.2018

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 1)

{
    print "Usage of GTH_rename_splitoutput.pl\n\n";
    print "perl GTH_rename_splitoutput.pl <gth_split_output>\n";
    print "where <gth_split_output> is the output file from running GenomeThreader with Nextflow splitting the queries\n";
    print "For example, >perl GTH_rename_splitoutput.pl gth_split_output.gff > gth_split_renamed.gff\n";
    exit;
}

my $gth_split_output = $ARGV[0];
my $line;
my @temp;
my $GeneID;

my $GeneChrom;
my $Genemethod;
my $Genefeature;
my $Genestart;
my $Geneend;
my $Genescore;
my $Genestrand;
my $Geneframe;
my $Genecomment;
my $RNAChrom;
my $RNAmethod;
my $RNAfeature;
my $RNAstart;
my $RNAend;
my $RNAscore;
my $RNAstrand;
my $RNAframe;
my $RNAcomment;
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
open(GTH,"$gth_split_output") || die "ERROR: cannot open $gth_split_output\n";

while (<GTH>) {
	$line = $_;
	chomp $line;
	
	@temp =split(/\t+/,$line);
	$feature = $temp[2];

	if ($feature eq "gene") {
	
		$GeneChrom = $temp[0];
		$Genemethod = $temp[1];
		$Genefeature = $temp[2];
		$Genestart = $temp[3];
		$Geneend = $temp[4];
		$Genescore = $temp[5];
		$Genestrand = $temp[6];
		$Geneframe = $temp[7];
		$Genecomment = $temp[8];
		
	} elsif ($feature eq "mRNA") {
		$RNAChrom = $temp[0];
		$RNAmethod = $temp[1];
		$RNAfeature = $temp[2];
		$RNAstart = $temp[3];
		$RNAend = $temp[4];
		$RNAscore = $temp[5];
		$RNAstrand = $temp[6];
		$RNAframe = $temp[7];
		$RNAcomment = $temp[8];
		
		($GeneID) =($RNAcomment =~/;Target=(\S+)\s/);
		print STDOUT $GeneChrom."\t".$Genemethod."\t".$Genefeature."\t".$Genestart."\t".$Geneend."\t".$Genescore."\t".$Genestrand."\t".$Geneframe."\tID=".$GeneID."\n";	
		print STDOUT $RNAChrom."\t".$RNAmethod."\tmatch_part\t".$RNAstart."\t".$RNAend."\t".$RNAscore."\t".$RNAstrand."\t".$RNAframe."\tID=".$GeneID."-RA;Parent=".$GeneID."\n";
	
	} else {
		$Chrom = $temp[0];
		$method = $temp[1];
		$feature = $temp[2];
		$start = $temp[3];
		$end = $temp[4];
		$score = $temp[5];
		$strand = $temp[6];
		$frame = $temp[7];
		$comment = $temp[8];
		print STDOUT $Chrom."\t".$method."\t".$feature."\t".$start."\t".$end."\t".$Genescore."\t".$strand."\t".$frame."\tID=".$GeneID."-RA\n";	
	}
}

close (GTH);

exit