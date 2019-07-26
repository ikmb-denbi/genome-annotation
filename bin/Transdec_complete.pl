#!/usr/bin/env perl

use strict;
use Getopt::Long;


my $usage = qq{
perl my_script.pl

Select only sequences with "type:complete" in the header.

  Getting help:
    [--help]
  Input:
    [--fasta_in filename]
		The name of the FASTA file to read. 
	[--gff_in filename]
		The name of the GFF file to read.
  Ouput:    
    [--fasta_out filename]
        The name of the output FASTA file. 
    [--gff_out filename]
    	The name of the output GFF file.
};

my $fasta_in = undef;
my $gff_in = undef;
my $fasta_out = undef;
my $gff_out = undef;
my $geneIDs_file = "Transdecoder_complete_GeneIDs.txt";
my $FASTA_OUT = undef;
my $IDS_FILE = undef;

my $header = undef;
my $type = undef;
my $hit = 0;
my $GeneID = undef;
my @completeIDs = undef;
my $counter = 0;
my $help;

my $help;

GetOptions(
    "help" => \$help,
    "fasta_in=s" => \$fasta_in,
    "gff_in=s" => \$gff_in,
    "fasta_out=s" => \$fasta_out,
    "gff_out=s" => \$gff_out);
    
    
# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($fasta_out) {
    open($FASTA_OUT, ">$fasta_out") or die("Cannot open $fasta_out");
}

if ($geneIDs_file) {
	open($IDS_FILE, ">$geneIDs_file") or die("Cannot open $geneIDs_file");
}

# open the Transdecoder output FASTA file 
open (my $FASTA_IN, '<', $fasta_in) or die "FATAL: Can't open file: $fasta_in for reading.\n$!\n";

while (<$FASTA_IN>) {

	my $line = $_;
	chomp $line;
	
	if ($line =~ /^>/) {
		($header) = $line;
		($type) = ($line=~/.*type:(\w+)\s.*/);
		($GeneID) = ($line=~/^>(\S+)\s.*/);
		if ($type eq "complete") {
			print $FASTA_OUT $line."\n";
			print $IDS_FILE $GeneID."\n";
			$hit = 1;
			$counter++;
			next;
		} else {
			$hit = 0;
			next;
		}
		
	} else {
		if ($hit == 1) {
			print $FASTA_OUT $line."\n";
			next;
		} else {
			next;
		}
	}
}

close $FASTA_IN;
close $FASTA_OUT;
close $IDS_FILE;

printf $counter." complete sequences.\n";

# 
my $cmd_grep = "grep -F -f $geneIDs_file $gff_in > $gff_out" ; 
system($cmd_grep);

