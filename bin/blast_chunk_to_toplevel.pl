#!/usr/bin/env perl
# A script to convert blast hits on a sub-partitioned assembly back to the top level coordinates

use strict;
use Getopt::Long;
use Data::Dumper;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--blast filename]
		The name of the BLAST report  to read. 
    [--agp filename]
		The name of AGP file to read

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $blast = undef;
my $agp = undef;
my $help;

GetOptions(
    "help" => \$help,
    "blast=s" => \$blast,
    "agp=s" => \$agp,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my @mappings;
my %scaffold_lengths;

open (my $AGP, '<', $agp) or die "FATAL: Can't open file: $agp for reading.\n$!\n";

foreach my $line ( <$AGP> ) {

	chomp($line);

	my ($chr,$start,$stop,$strand,$phase,$chunk_name,$chunk_start,$chunk_stop,$chunk_strand) = split("\t", $line);
	
	push @mappings, { "chr" => $chr, "start" => $start, "stop" => $stop, "phase" => $phase, "chunk_name" => $chunk_name, "chunk_start" => $chunk_start,"chunk_stop" => $chunk_stop, "chunk_strand" => $chunk_strand };

	$scaffold_lengths{$chr} = $stop;
}

close($AGP);

open (my $BLAST, '<', $blast) or die "FATAL: Can't open file: $blast for reading.\n$!\n";

my $agp_info = undef;

foreach my $line ( <$BLAST> ) {

	chomp($line);
	# qseqid sseqid sstart send slen pident qlen qstart qend length mismatch gapopen evalue bitscore
	my ($qseqid,$sseqid,$sstart,$send,$slen,$pident,$qlen,$qstart,$qend,$length,$mismatch,$gapopen,$evalue,$bitscore) = split("\t", $line);

	if (!defined $agp_info || $agp_info->{'chunk_name'} ne $qseqid) {
		($agp_info) = grep { $qseqid eq $_->{chunk_name} } @mappings;
	}

	#printf Dumper(\$agp_info) . "\n";

	die "Ooops, I couldn't find the AGP mapping for $qseqid...\n" if !$agp_info ;
	#printf "Offset for $qseqid is " . $agp_info->{"start"} . "\n";	

	my $actual_start = $agp_info->{"start"}+$qstart ;
	my $actual_stop = $agp_info->{"start"}+$qend ;
	my $chr_length = $scaffold_lengths{$agp_info->{"chr"}} ;
	printf $agp_info->{"chr"} . "\t" . $sseqid . "\t" . $sstart . "\t" . $send . "\t" . $slen . "\t" . $pident . "\t" . $chr_length ."\t" . $actual_start . "\t" . $actual_stop . "\t". $length . "\t" . $mismatch . "\t" . $gapopen . "\t" . $evalue . "\t". $bitscore . "\n" ;
}

close($BLAST);

exit 0;
