#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the GFF file to read. 

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gff = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my @transcripts;

open (my $GFF, '<', $gff) or die "FATAL: Can't open file: $gff for reading.\n$!\n";

while (<$GFF>) {

	chomp;
	my $line = $_;

	my $info = (split "\t",$line)[-1];
	my %data;
	foreach my $element (split ";",$info) {
		my ($key,$value) = split("=",$element);
		if ($key eq "Target") {
			my $transcript_id = (split " ",$value)[0];
			printf $transcript_id . "\n";
		}
	}
}

close($GFF);

