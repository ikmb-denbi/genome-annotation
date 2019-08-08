#!/usr/bin/env perl
# This file writes a valid pasa config file from a template since we need the actual full path to the work 
# sub directory for nextflow to properly track the sqlite file. This can perhaps be replaced with some awk magic.

use strict;
use Getopt::Long;
use Cwd;

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
my $wd = getcwd ;

foreach my $line (<$IN>) {

	chomp($line);
	if ($line =~ /^DATABASE.*/) {
		printf "DATABASE=${wd}/pasa_DB.sqlite\n";
	} else {
		printf $line . "\n" ;
	}
}

close($IN);
