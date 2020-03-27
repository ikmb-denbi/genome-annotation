#!/usr/bin/env perl
# A script to convert a BLAST report to a target list for Exonerate

use strict;
use Getopt::Long;
use Data::Dumper;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
	[--max_intron_size]
		Maximum length of intron for spliced alignments
	[--min_bit]
		Minimum bitscore for match to be considered valid (default: 50)
	[--min_id]
		Miniumum id% for match to be considered valid (default: 0.6)
	[--length_percent ]
		Minimum length percentage for the match to be considered valid (default: 0.6)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my %targets;
my $min_bit = 25.0; # bit score of the match must be at least this high
my $min_id = 50.0; # match must be at least this similar
my $length_percent = 0.6; # At least this much of the query has to align
my $max_intron_size = undef;

my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "max_intron_size=i" => \$max_intron_size,
    "min_bit=i" => \$min_bit,
    "min_id=f" => \$min_id,
    "length_percent=i" => \$length_percent,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (!defined $max_intron_size){
	die "Must provide a maximum intron size (--max_intron_size)\n" ;
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

my %bucket;

while(<$IN>) {
	
	chomp; 
	my $line = $_; 
	# I       Y74C9A.5.1      120     344     474     97.778  50000   30443   79768   225     5       0       2.72e-146       467
	# chr1    sp|Q14393|GAS6_HUMAN    278     386     678     61.4    191232075       5657    5184    158     12      1       2.1e-46 189.5

	my ($qseqid, $sseqid, $sstart, $send, $slen, $pident, $qlen, $qstart, $qend, $length, $mismatch, $gapopen, $evalue, $bitscore) = split("\t", $line);

	my %entry = (
		"query_id" => $qseqid,
		"target_id" => $sseqid,
		"target_start" => $sstart,
		"target_end" => $send,
		"target_length" => $slen,
		"query_length" => $qlen,
		"query_start" => $qstart,
		"query_end" => $qend,
		"aln_length" => $length,
		"pident" => $pident,
		"mismatch" => $mismatch,
		"gapopen" => $gapopen,
		"evalue" => $evalue,
		"bitscore" => $bitscore
	);
	# If this match is in reverse direction
	if ($entry{"target_end"} < $entry{"target_start"} ) {
		my $new_start = $entry{"target_end"};
		my $new_end = $entry{"target_start"};
		$entry{"target_start"} = $new_start;
		$entry{"target_end"} = $new_end;
	}
	if ($entry{"query_end"} < $entry{"query_start"}){
		 my $new_start = $entry{"query_end"};
                my $new_end = $entry{"query_start"};
                $entry{"query_start"} = $new_start;
                $entry{"query_end"} = $new_end;

	}

	# add to our inventory ( a hash of hashes using the query and target id as keys)
	if ( exists $bucket{$sseqid} ) {
		if (exists $bucket{$sseqid}{$qseqid}) {
			push( @{ $bucket{$sseqid}{$qseqid} }, \%entry ); 
		} else {
			$bucket{$sseqid}{$qseqid} = [ \%entry ];
		}
	} else {
		$bucket{$sseqid}{$qseqid} = [ \%entry ] ;
	}

}

# All BLAST entries are grouped by query and target sequence, now stitch into clusters based on max_intron length

# Iterate over each protein to build clusters
foreach my $query ( keys %bucket ) {

	#printf "Checking matches for $query\n";	

	my $data = $bucket{$query};

	# every protein - chromosome/scaffold relationship
	# Our target are the outer most mapping coordinates +/- twice the max intron size
	# Iterating over each mapped scaffold...
	foreach my $target ( keys %$data ) {

		#printf "Building data for target $target\n";
		
		# Sort the query/target specific matches based on their starting position
		my $matches = $bucket{$query}{$target};

		# Coordinates along the genomic scaffolds
		my @sorted_matches = sort { $a->{"query_start"} <=> $b->{"query_start"} } @{$matches};
		# Coordinates along the target protein
		my @sorted_query_matches = sort { $a->{"target_start"} <=> $b->{"target_start"} } @{$matches};


		# Iterate over sorted genomic intervals; merge into one region for exonerate
		# break into multiple regions if any given gap is larger than 4* the max intron length

		my @clean_matches;
		foreach my $match (@sorted_matches) {

			if ( $match->{'bitscore'} >= $min_bit ) {
				push(@clean_matches,$match);
			} else {
				my $match_name = $match->{'target_id'} . ":" . $match->{'query_id'} . ":" . $match->{'query_start'} . "-" . $match->{'query_end'} . ":" . $match->{'bitscore'};
				print STDERR "Skipping match $match_name due to low score...\n";
			}
		}

		my $first_entry = shift @clean_matches;

		my $this_start = $first_entry->{"query_start"};
		my $this_end = $first_entry->{"query_end"};

		my @bucket;
		push(@bucket,$first_entry);

		# stitch matches into clusters, gaps may be no longer than 2 times the max intron size
		# This prevents creating absurdely large clusters in case a proteins maps to a scaffold more than once
		foreach my $match (@clean_matches) {
			if ($match->{"query_start"} > $this_start+($max_intron_size*2) ) {
				# this gap is too large - dump out what we have and start new cluster
				print_cluster(\@bucket,$max_intron_size);
				@bucket = ();
			} else {
				push(@bucket,$match);
			}
	
		}
		print_cluster(\@bucket,$max_intron_size);

	}		

}

sub print_cluster {

	my @matches = @{$_[0]};
	my $max_intron_size = $_[1];

	if (scalar @matches > 0) {

		#print STDERR "Building evidence cluster with intron size $max_intron_size\n";

		my $first_entry = @matches[0];
		my $last_entry = @matches[-1];
	
		my $this_start = $first_entry->{'query_start'} ;
		my $this_end = $last_entry->{'query_end'};
	
		# Define the final genomic coordinates adding twice the maximum intron size as flanking regions (or start/end of scaffold, whichever comes first)
        	my $target_start = ($this_start <=  $max_intron_size*2) ? 1 : ($this_start-$max_intron_size*2);
	        my $target_end = ($this_end+$max_intron_size*2) >= $first_entry->{"query_length"} ? $first_entry->{"query_length"} : ($this_end+$max_intron_size*2);

        	printf $first_entry->{"target_id"} . "\t" . $first_entry->{"query_id"} . "\t" . $target_start . "\t" . $target_end . "\n";
	}

}
