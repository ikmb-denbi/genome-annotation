#!/usr/bin/env perl

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
my $min_id = 60; # match must be at least this similar
my $length_percent = 0.6; # At least this much of the query has to align
my $max_intron_size = undef;

my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "max_intron_size=i" => \$max_intron_size,
    "min_bit=i" => \$min_bit,
    "min_id=i" => \$min_id,
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
	if ( exists $bucket{$qseqid} ) {
		if (exists $bucket{$qseqid}{$sseqid}) {
			push( @{ $bucket{$qseqid}{$sseqid} }, \%entry ); 
		} else {
			$bucket{$qseqid}{$sseqid} = [ \%entry ];
		}
	} else {
		$bucket{$qseqid}{$sseqid} = [ \%entry ] ;
	}

}

# All BLAST entries are grouped by query and target sequence, now stitch into clusters based on max_intron length

# Iterate over each protein to build clusters
foreach my $query ( keys %bucket ) {

	my $data = $bucket{$query};

	# every protein - chromosome/scaffold relationship
	# Our target are the outer most mapping coordinates +/- twice the max intron size
	foreach my $target ( keys %$data ) {
		
	
		# Sort the query/target specific matches based on their starting position
		my $matches = $bucket{$query}{$target};

		my @sorted_matches = sort { $a->{"target_start"} <=> $b->{"target_start"} } @{$matches};
		my @sorted_query_matches = sort { $a->{"query_start"} <=> $b->{"query_start"} } @{$matches};

		# Sorted by query (protein) positions to determine how much of the query sequence is aligned in this scaffold
		my $first_query_entry = @sorted_query_matches[0];
		my $last_query_entry = @sorted_query_matches[-1];

		my $query_length = $first_query_entry->{"query_length"};
		my $match_length = $last_query_entry->{"query_end"} - $first_query_entry->{"query_start"};
		my $fraction = $match_length/$query_length;
	
		# Sorted by target positions to determine the genomic boundaries for subsequent exonerate alignments
		my $first_entry = @sorted_matches[0];
		my $last_entry = @sorted_matches[-1];

		# is this overall a valid match?
		if ($first_entry->{"bitscore"} >= $min_bit && $last_entry->{"bitscore"} >= $min_bit && $fraction >= $length_percent ) {

			my $this_start = $first_entry->{"target_start"} ;
			my $this_end = $last_entry->{"target_end"} ;

			my $target_start;
			my $target_end; 

			# Define the final genomic coordinates adding twice the maximum intron size as flanking regions (or start/end of scaffold, whichever comes first)
			my $target_start = ($this_start <=  $max_intron_size*2) ? 1 : ($this_start-$max_intron_size*2);
                	my $target_end = ($this_end+$max_intron_size*2) >= $first_entry->{"target_length"} ? $first_entry->{"target_length"} : ($this_end+$max_intron_size*2);
		
			 printf $first_entry->{"query_id"} . "\t" . $first_entry->{"target_id"} . "\t" . $target_start . "\t" . $target_end . "\n";
		}

	}

}

