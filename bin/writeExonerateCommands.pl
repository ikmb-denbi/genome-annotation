#!/bin/wnv perl

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
my $source = undef;
my %source_keys = (
	"proteins" => "protein2genome",
	"transcripts" => "est2genome"	
)

my $help;

GetOptions(
    "help" => \$help,
    "matches=s" => \$matches,
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

open (my $QUERIES, '<', $queries) or die "FATAL: Can't open file: $queries for reading.\n$!\n";

my $analysis = $source_keys{$source};

while (<$QUERIES) {

	$line = $_;
	chomp $line;

	@temp = split(/\t+/,$line);
	$query = $temp[0];
    	$target = $temp[1];
	
	my $query_clean = $query;
	$query_clean =~ s/\|/_/g ;

	my $fa_clean = $query_clean + ".fa" ;

	printf "exonerate --model $analysis --softmasktarget yes --bestn 1 --minintron 20 --maxintron 20000  --showalignment false --showtargetgff true $fa_clean $genome_basename/$target.fa > $fa_clean.exonerate.out\n";
}


