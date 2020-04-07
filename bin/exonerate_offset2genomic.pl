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

print STDERR "Starting Exonerate offset2genomic parser...\n";

printf "## from $infile\n";

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

while (<$IN>) {
	
	chomp; 
	my $line = $_; 
	
	if ($line =~ ".*exonerate\:.*") {
		my ($target,$method,$feature,$start,$stop,$score,$strand,$frame,$attributes) = split("\t",$line);
		my @info = split(/[:,-]/,$target);
		my $offset = @info[1] ;
		my $target_base = @info[0];
		my $adjusted_start = $start+$offset-1;
		my $adjusted_stop =$stop+$offset-1 ;
		my $method = "protein2genome" ;
		printf "$target_base\t$method\t$feature\t$adjusted_start\t$adjusted_stop\t$score\t$strand\t$frame\t$attributes\n";
	}
}

close $IN;

print STDERR "Finished converting back to genomic coordinates\n";
