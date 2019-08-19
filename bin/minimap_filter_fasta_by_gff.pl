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
    [--fasta filename]
		The name of the fasta file

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gff = undef;
my $fasta = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "fasta=s" => \$fasta,
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
			unless (grep( /^$transcript_id$/, @transcripts ) ) {
				push(@transcripts,$transcript_id);
			}
		}
	}
}

close($GFF);

open (my $FASTA, '<', $fasta) or die "FATAL: Can't open file: $fasta for reading.\n$!\n";

my $skip = 1;

while (<$FASTA>) {

	chomp;
	my $line = $_;

	if ($line =~ /^>.*/) {
		my $definition = (split ">",$line)[-1];
		my $transcript_id = (split " ",$definition)[0];
		if ( grep( /^$transcript_id$/, @transcripts ) ) {
			$skip = 0;
		} else {
			$skip = 1;
		}
	}

	if ($skip == 0) {
		printf $line ."\n";
	}
		
}

close ($FASTA);
