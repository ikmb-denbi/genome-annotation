#!/usr/bin/perl

use strict;
use Getopt::Long;
use POSIX;

my $usage = qq{
perl split_training.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
	The name of the input file (in Genbank format) with peptide sequences to train Augustus.
    [--percent integer]
	Percentage from total peptides to go to training (default=90)
};

my $infile =undef;
my $percent = undef;
my $count = 0;
my $trainint = undef;
my $testint = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "percent=i" => \$percent);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (!defined $percent){
        $percent = 90;
} else {
	unless (($percent <= 100) && ($percent > 0)) {
		die "Wrong percentage value\n";
	}
}
 
open(my $IN, '<', $infile) or die "FATAL: Cannot open file: $infile for reading.\n$!\n";

while (<$IN>) {
	
	chomp; 
	my $line = $_; 
	if ($count == 0) {
		unless ($line =~ /LOCUS/) {
			die "Input file is not in GenBank format.\n";
		}
	}
	
	if ($line =~ /LOCUS/) {
		$count ++;
	}	
}

close($IN);

# How many sequences to training?
$trainint = ceil($count * ($percent/100));
$testint = $count - $trainint;

# Get random splitting of gene set
my $cmd_run = "randomSplit.pl $infile $testint\n";

print "Running randomSplit.pl on ".$infile.".\nTraining set: ".$trainint." peptides.\nTest set: ".$testint.".\n";

system($cmd_run);
