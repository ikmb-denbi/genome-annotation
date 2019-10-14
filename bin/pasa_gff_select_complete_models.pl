#!/usr/bin/env perl
Return transcript models flagged as complete by PASA

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

my $valid = 0;

while (<$IN>) {

	chomp;
	my $line = $_;

	if ($line =~ /^#.*/) {
		next;
	}
	my @elements = split "\t", $line;

	my %entry = ( "seq_name" => $elements[0], "source" => $elements[1], "feature" => $elements[2], "start" => $elements[3], "stop" => $elements[4],
		"score" => $elements[5], "strand" => $elements[6], "phase" => $elements[7], "attributes" => $elements[8]);

	if ( $entry{'feature'} eq "gene" ) {
		if ($entry{'attributes'} =~ /.*complete.*/ ) {
			$valid = 1;
		} else {
			$valid = 0;
		}
	}

	if ($valid == 1) {
		printf $line . "\n";
	}
}

close($IN);
