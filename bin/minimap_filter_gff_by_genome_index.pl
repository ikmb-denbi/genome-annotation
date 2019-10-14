#!/usr/bin/env perl
# Return Minimap alignments overlapping regions defined in a fasta index
use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the GFF file to read. 
    [--index filename]
		The name of the genome fasta index .fai

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gff = undef;
my $index = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "index=s" => \$index,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $INDEX, '<', $index) or die "FATAL: Can't open file: $index for reading.\n$!\n";

my @scaffolds ;

while (<$INDEX>) {

	chomp;
	my $line = $_;

	my $seq_name = (split "\t",$line)[0];

	push(@scaffolds,$seq_name);
}

close($INDEX);

open (my $GFF, '<', $gff) or die "FATAL: Can't open file: $gff for reading.\n$!\n";

while (<$GFF>) {

	chomp;
	my $line = $_;

	my $seq_name = (split "\t",$line)[0];
	if ( grep( /^$seq_name$/, @scaffolds )  ) {
		printf $line ."\n";
	}
}

close($GFF);
