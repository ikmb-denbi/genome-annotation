#!/usr/bin/env perl
# GFF from Minimap to hints format for augustus

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
my $source = "minimap";
my $GeneID = undef;
my $help;

my $src = "T";
my $pri = 4;
my $hintfeature = "exonpart";

my $help;

GetOptions(
    "help" => \$help,
    "source=s" => \$source,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# open the minimap GFF file
open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

while (<$IN>) {
	
	my $line = $_;
	chomp $line;

	if ( $line =~ /^#.*/ ) {
		printf $line . "\n";
	} else {
	#I       genome  cDNA_match      15053797        15053899        100.0   -       .       ID=F31C3.6a.p1;Target=F31C3.6a 1560 1662
		my ($chr,$origin,$feature,$from,$to,$score,$strand,$phase,$info) = split("\t",$line);
		my $group = (split /[;,=]/ , $info)[1];
		printf $chr . "\t" . $source . "\t" . $hintfeature . "\t" . $from . "\t" . $to . "\t" . $score . "\t" . $strand . "\t" . $phase . "\t" . "src=$src;grp=$group;pri=$pri\n"

	}
}

close $IN;
