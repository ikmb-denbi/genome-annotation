#!/usr/bin/env perl
# Run exonerate based on regions detected with blast

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--matches filename]
		The name of the match file to read. 
	[--assembly_index filename]
		Name of the CDBtools genome assembly index
	[--query_index filename]
		Name of the CDBtools protein fasta index
	[--analysis name]
		Exonerate algorithm to use (est2genome or protein2genome)
	[--max_intron_size]
		Maximum length of intron to consider for spliced alignments
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $matches = undef;
my $assembly_index = undef;
my $max_intron_size = undef;
my $query_index = undef;
my $analysis = undef;
my $help;

GetOptions(
    "help" => \$help,
    "matches=s" => \$matches,
    "assembly_index=s" => \$assembly_index,
    "max_intron_size=i" => \$max_intron_size,
    "query_index=s" => \$query_index,
    "analysis=s" => \$analysis,
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

print STDERR "Preparing exonerate jobs!\n" ;

open (my $IN, '<', $matches) or die "FATAL: Can't open file: $matches for reading.\n$!\n";

while (<$IN>) {
	
	chomp; 
	my $line = $_; 
	
	# A pre-filtered BLAST match	
	# get all relevant columns
	my ($qseqid,$sseqid,$target_start,$target_end) = split("\t", $line);
	
	# get query sequence
	my $query_clean = $qseqid;
	$query_clean =~ s/\|/_/g ;
	my $fa_clean = "$query_clean._query_.fa" ;
	
	# THIS SHOULD COME FROM ANOTHER SCRIPT
	#my $cmd_query = "cdbyank $query_index -a '$qseqid' > $fa_clean" ;
	#system($cmd_query);
	
	# get subregion of target scaffold
	my $region = $sseqid . ":" . $target_start . "-" . $target_end ;
	my $region_name = $sseqid . "_" . $target_start . "_" . $target_end ;
	my $cmd_target = "samtools faidx $assembly_index $region > $region_name._target_.fa" ;
	system($cmd_target);
	
	# Run exonerate on these data
	my $cmd_run = "exonerate --model $analysis --softmasktarget --percent 20 --bestn 1 --minintron 20 --maxintron $max_intron_size  --showalignment false --showtargetgff true $fa_clean $region_name._target_.fa > subjob_$query_clean.$region_name.exonerate.align\n";
	
	printf($cmd_run);

}

close($IN);

print STDERR "Finished building exonerate jobs!\n";

