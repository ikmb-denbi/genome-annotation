#!/usr/bin/env perl
# Converts Exonerate output to GFF3 hints
use strict;
use Getopt::Long;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
    [--source string]
		A valid source for processing (est, protein or trinity)
    [--pri int]
		Priority for protien-derived hints (default: 5)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $source = "protein2genome";
my $GeneID = undef;
my $pri = 5;
my $src = "P";
my $help;


GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "source=s" => \$source,
    "pri=i", => \$pri,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}


# open the exonerate report 
open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

while (<$IN>) {
	
	my $line = $_;
	chomp $line;

	my ($Chrom,$met,$feature,$start,$end,$score,$strand,$frame,$comment) = split(/\t+/,$line);

	if ($feature eq "gene") {
		($GeneID) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)\s;\s/);
	} elsif ($feature eq "exon") {
		printf $Chrom."\t".$source."\texonpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} elsif ($feature eq "cds") {
		printf $Chrom."\t".$source."\tCDSpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} elsif ($feature eq "intron") {
		printf $Chrom."\t".$source."\tintronpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} elsif ($feature eq "utr5" || $feature eq "utr3") {
		printf $Chrom."\t".$source."\tUTRpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} else {
		next;
	}
}

close $IN;
