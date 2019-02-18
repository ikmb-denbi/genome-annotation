#!/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--matches filename]
		The name of the tabular match file to read. 
    [--target_root ]
		The root directory where the targets live
    [--source]
		The type of query (protein,transcript)
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $matches = undef;
my $db = undef;
my $target_root = undef;
my $max_intron_size = 20000;
my $source = undef;
my %source_keys = ( "proteins" => "protein2genome","transcripts" => "transcript2genome", "est" => "est2genome" );

my $help;

GetOptions(
    "help" => \$help,
    "matches=s" => \$matches,
    "source=s" => \$source,
    "max_intron_size=i" => \$max_intron_size,
    "db=s" => \$db,
    "target_root=s" => \$target_root,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

if (!$source) {
	exit 1; "Must provide a source name (proteins or transcripts)\n" ;
}

open (my $QUERIES, '<', $matches) or die "FATAL: Can't open file: $matches for reading.\n$!\n";

my $analysis = $source_keys{$source};

while (<$QUERIES>) {

	my $line = $_;
	chomp $line;

	my @temp = split(/\t+/,$line);
	my $query = $temp[0];
    	my $target = $temp[1];
	
	my $query_clean = $query;
	$query_clean =~ s/\|/_/g ;

	my $fa_clean = "$query_clean.fa" ;

	printf "exonerate --model $analysis --softmasktarget yes --bestn 1 --minintron 20 --maxintron $max_intron_size  --showalignment false --showtargetgff true $fa_clean $target_root/$target.fa > $fa_clean.$target.exonerate.out\n";
}


