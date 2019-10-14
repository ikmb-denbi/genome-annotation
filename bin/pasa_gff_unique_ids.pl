#!/usr/bin/env perl
# Merge PASA cDNA alignments across chunks and make new valid IDs

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;

my $help;

GetOptions(
    "help" => \$help,
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

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";


while (<$IN>) {

	my $line = $_;
	chomp $line;

	if ($line =~ /^#.*/) {

		printf $line . "\n";

	} else {
		my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

		my @attributes ;

		foreach my $attrib (split ";",$info) {

			my ($key,$value) = split("=",$attrib);
			if ($key eq "ID" || $key eq "Parent") {
				$value .= "-$seq" ;
			}
			push(@attributes,"$key=$value");
		}

		printf $seq . "\t" . $source . "\t" . $feature . "\t" . $start . "\t" . $stop . "\t" . $phase . "\t" . $strand . "\t" . $score . "\t" . join(";",@attributes) . "\n";

	}
}
