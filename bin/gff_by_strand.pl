#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
    [--max_intron int]
		The maximum size of intron to consider (default: 20000)

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $max_intron = 20000;

my $help;

GetOptions(
    "help" => \$help,
    "max_intron=i" => \$max_intron,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}


open (my $PLUS, '>', "hints.plus.gff");
open (my $MINUS, '>', "hints.minus.gff");

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";


while (<$IN>) {

	my $line = $_;
	chomp $line;

	my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

	if ($strand eq "+") {
		printf $PLUS $line . "\n";
	} elsif ($strand eq "-") {
		printf $MINUS $line . "\n";
	}
}

close($IN);
close($PLUS);
close($MINUS);


