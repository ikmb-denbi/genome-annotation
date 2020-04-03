#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--genome_fai filename]
		The name of the assembly fasta index
    [--bed filename]
		A bed file of target regions
    [--hints filename]
		A hints file in GFF format
    [--aug_conf filename]
		An augustus custom config file with hint weights
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

# augustus_from_regions.pl --genome_fai $genome_fai --bed $regions --hints $hints --aug_conf $AUG_CONF --isof $params.isof --utr $params.UTR
my $outfile = undef;
my $infile = undef;
my $genome_fai = undef;
my $bed = undef;
my $hints = undef;
my $aug_conf = undef;
my $isof = undef;
my $utr = undef;
my %dictionary;
my $model = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "model=s" => \$model,
    "genome_fai=s" => \$genome_fai,
    "hints=s" => \$hints,
    "utr=s" => \$utr,
    "isof=s" => \$isof,
    "bed=s" => \$bed,
    "aug_conf=s" => \$aug_conf,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# Read FAI file to get list of scaffolds

open (my $FAI, '<', $genome_fai) or die "FATAL: Can't open file: $genome_fai for reading.\n$!\n";

while (<$FAI>) {

	my $line = $_;
	chomp($line);

	my ($chr,$len,$a,$b,$c) = split("\t",$line);

	$dictionary{$chr} = $len;
	
}

close($FAI);

# Go over the regions from the bed file and define jobs

open (my $BED, '<', $bed) or die "FATAL: Can't open file: $bed for reading.\n$!\n";

while (<$BED>) {

	my $line = $_;
	chomp($line);

	my ($chr,$from,$to,$strand) = split("\t",$line);

	# If this scaffold is part of our target chunk
	if ( exists $dictionary{$chr}) {
		my $scaffold_len = $dictionary{$chr};

		my $max_len = ($to > $scaffold_len) ? $to = $scaffold_len : $to = $to;

		my $outfile = $chr . "_" . $from . "-" . $to . "." . $strand . ".augustus.gff" ;
		my $infile = $chr . ".fa" ;

		my $direction = "";
		if ($strand eq "+") {
			$direction = "--strand=forward";
		} elsif ($strand eq "-") {
			$direction = "--strand=backward";
		}
		
		my $command = "augustus --species=$model $direction --alternatives-from-sampling=false --alternatives-from-evidence=false --hintsfile=$hints --gff3=on --UTR=$utr --alternatives-from-evidence=$isof --extrinsicCfgFile=$aug_conf --hintsfile=$hints --predictionStart=$from --predictionEnd=$to --uniqueGeneId=true $chr.fa > $outfile" ;
		printf $command . "\n" ;
	}
}
close($BED);

