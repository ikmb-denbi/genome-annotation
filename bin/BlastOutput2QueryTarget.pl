#!/usr/bin/env perl

# By M. Torres on 23.10.2017

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 3)
{
    print "Usage of BlastOutput2QueryTarget.pl\n\n";
    print "perl BlastOutput2Sequences.pl <blast_output> <eval_cutoff> <output_file>\n";
    print "where <blast_output> is the input file with blast results,\n";
    print "      <eval_cutoff> is the evalue cutoff to use for blast matches,\n";
    print "      <output_file> is the output file with the query-target combinations,\n";
    print "For example, >perl BlastOutput2Sequences.pl blast_output.txt 0.05 blast_query_subject.txt\n";
    exit;
}


my $blast_output		= $ARGV[0]; # BLAST OUTPUT FILE
my $evalue_cutoff		= $ARGV[1]; # EVALUE CUTOFF
my $output_file			= $ARGV[2]; # OUTPUT FILE
my $line;
my @temp;
my $query;						# QUERY PROTEIN OF THAT HIT
my $subject;					# TARGET SCAFFOLD OF THAT HIT
my $evalue;						# EVALUE OF THAT HIT
   

# OPEN FILE WITH BLAST RESULTS:   
open(BLAST,"$blast_output") || die "ERROR: cannot open $blast_output\n";

# OPEN FILE TO WRITE OUTPUT ON:
open(OUTFILE,">$output_file") || die "ERROR: cannot open $output_file\n";

# ITERATE OVER BLAST RESULTS FILE:
while(<BLAST>)
{
	$line                  = $_;
	chomp $line;
	
	# IGNORE COMMENT LINES IF ANY:
	if (substr($line,0,1) ne '#')
	{
    	@temp                  = split(/\t+/,$line);
        # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        $query				 = $temp[0];
        $subject            = $temp[1];
        $evalue             = $temp[10];
        
        # ONLY WHEN EVAL IS EQUAL OR SMALLER THAN THE CUTOFF:
       	if ($evalue <= $evalue_cutoff)
       	{
           	# PRINT THE QUERY AND THE TARGET SCAFFOLD/CROMOSOME ON THE OUTPUT FILE:
           	print OUTFILE "$query\t$subject\n";
        }
	}
}
close(BLAST);
close(OUTFILE);