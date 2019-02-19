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
	[--max_intron_size]
		Maximum length of intron for spliced alignments
	[--min_bit]
		Minimum bitscore for match to be considered valid (default: 50)
	[--min_id]
		Miniumum id% for match to be considered valid (default: 0.6)
	[--length_percent ]
		Minimum length percentage for the match to be considered valid (default: 0.6)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my %targets;
my $min_bit = 50.0; # bit score of the match must be at least this high
my $min_id = 60; # match must be at least this similar
my $length_percent = 0.6; # At least this much of the query has to align
my $max_intron_size = undef;

my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
	"max_intron_size=i" => \$max_intron_size,
	"min_bit=i" => \$min_bit,
	"min_id=i" => \$min_id,
	"length_percent=i" => \$length_percent,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (!defined $max_intron_size){
	die "Must provide a maximum intron size (--max_intron_size)\n" ;
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

while(<$IN>) {
	
	chomp; 
	my $line = $_; 
	my ($qseqid, $sseqid, $sstart, $send, $slen, $pident, $qlen, $length, $mismatch, $gapopen, $evalue, $bitscore) = split("\t", $line);

	my $target_start;
	my $target_end;
	
	# Get the alignment region +/- twice the maximium intron length
	if ($length >= $qlen*$length_percent && $bitscore >= $min_bit && $pident >= $min_id) {
		$target_start = ($sstart <=  $max_intron_size*2) ? 1 : ($sstart-$max_intron_size*2);
		$target_end = ($send+$max_intron_size*2) >= $slen ? $slen : ($send+$max_intron_size*2);
		printf "$qseqid\t$sseqid\t$target_start\t$target_end\n";
	}
		
}
