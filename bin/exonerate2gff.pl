#!/bin/env perl

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
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $source = undef;
my $GeneID = undef;
my $help;

my %source_keys = (
	"est" => {
		"src" => "T",
		"pri" => 4,
		"source" => "est2genome"
	},
	"protein" => {
		"src" => "P",
		"pri" => 3,
		"source" => "protein2genome"
		
	},
	"trinity" => {
		"src" => "T",
		"pri" => 4,
		"source" => "trinity2genome"
	}
) ;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "source=s" => \$source,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# make sure the source is valid
my @sources = keys %source_keys;

unless ( grep( /^$source$/, @sources ) ) {
  exit 1, "Did not provide a valid source (${join(',', @{ keys %source_keys })})\n";
}

# open the exonerate report 
open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

my $src = $source_keys{$source}{"src"};
my $pri = $source_keys{$source}{"pri"};
my $method = $source_keys{$source}{"source"};

while (<$IN>) {
	
	my $line = $_;
	chomp $line;

	my ($Chrom,$met,$feature,$start,$end,$score,$strand,$frame,$comment) = split(/\t+/,$line);

	if ($feature eq "gene") {
		($GeneID) =($comment =~/gene_id\s\w+\s;\ssequence\s(\S+)\s;\s/);
	} elsif ($feature eq "exon") {
		printf $Chrom."\t".$method."\texonpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;src=$pri\n";
	} elsif ($feature eq "cds") {
		printf $Chrom."\t".$method."\tCDSpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;pri=$pri\n";
	} elsif ($feature eq "intron") {
		printf $Chrom."\t".$method."\tintronpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;src=$pri\n";
	} elsif ($feature eq "utr5" || $feature eq "utr3") {
		printf $Chrom."\t".$method."\tUTRpart\t".$start."\t".$end."\t".$score."\t".$strand."\t".$frame."\tgrp=".$GeneID.";src=$src;src=$pri\n";
	} else {
		next;
	}
}

close $IN;
