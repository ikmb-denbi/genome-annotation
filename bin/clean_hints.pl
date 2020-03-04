#!/usr/bin/env perl

# Remove single-exon EST hints and low-coverage intron hints from rnaseq

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
my $min_cov = 10;
my $previous_id = undef;

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

my @bucket;

while (<$IN>) {

	my $line = $_;
	chomp $line;

	my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

	my %attribs;
	my @data = (split ";", $info) ;
	foreach my $a (@data) {
		my ($key,$value) = split("=",$a);
		$attribs{$key} = $value;
	}

	my $this_id = $attribs{'grp'} ;

	# Intron hints are only valid if they have a coverage equal or higher to 10X
	if ($source eq 'b2h') {
		if (defined $attribs{'mult'}) {
			my $multi = $attribs{'mult'};
			if ($multi >= $min_cov) {
				print $line . "\n";
			}
		}

	# EST hints are removed if they are single-exon (unspliced)
	} elsif ($source eq 'est2genome' || $source eq 'minimap_est' ) {
		if (defined $previous_id) {

			# A spliced hint continues
			if ($previous_id eq $this_id) {
				push(@bucket,$line);

			# A new ID starts, this hint group is complete 
			# Now check if multi-exon or not
			} else {
				if (scalar @bucket > 1) {
					foreach my $valid (@bucket) {
						print $valid . "\n";
					}
				}
				@bucket = ( $line );
			}
			# Set the new running ID
			$previous_id = $this_id ;

		} else {

			$previous_id = $this_id ;
			push(@bucket,$line);
		}

	# Proteins are fine, we don't filter those
	} else {
		printf $line . "\n";
	}

}

# empty bucket if not empty and more than one element
if (scalar @bucket > 1) {

	foreach my $valid (@bucket) {
        	print $valid . "\n";
        }

}
