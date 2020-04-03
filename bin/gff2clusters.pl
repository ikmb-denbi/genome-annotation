#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
    [--max_intron int]
		The maximum size of intron to consider (default: 20000)

  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $max_intron = 20000;

my $help;

GetOptions(
    "help" => \$help,
    "max_intron=i" => \$max_intron,
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

my $this_start = 0;
my $previous_start = 0;
my $previous_strand = "-";

my %this_block = undef;

while (<$IN>) {

	my $line = $_;
	chomp $line;

	my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

	# printf $seq . "\t" . $start .  "\n";

	if ($strand eq ".") {
		$strand = $previous_strand;
	}

	$this_start = $start;

	if (!defined $this_block{'seq_name'}) {
		#printf "Starting a new block\n";
                %this_block = ( 'seq_name' => $seq, 'start' => $start, 'stop' => $stop , 'strand' => $strand);

	# Existing block, same sequence - extend or finish
	} elsif ( $this_block{'seq_name'} eq $seq ) {
		
	        die "File not coordinate sorted! Previous start was $previous_start and this start is $this_start\n" if ($this_start < $previous_start );
	
		# Within range of the previous range and on the same strand
		if ($this_block{'stop'} >= ($start-$max_intron) && $strand eq $previous_strand) {
			$this_block{'stop'} = $stop;
		# out of range
		} else {

			#my $gap = ($start - $this_block{'stop'});
			#printf "Starting a new block - gap of $gap\n";

			# finish and print block
			$this_block{'start'} = 1 if ($this_block{'start'} - $max_intron < 0);
			$this_block{'stop'} = ($this_block{'stop'}+$max_intron);
			printf $this_block{'seq_name'} . "\t" . $this_block{'start'} . "\t" . $this_block{'stop'} . "\t" . $this_block{'strand'} . "\n";

			# start a new block
	                %this_block = ( 'seq_name' => $seq, 'start' => $start, 'stop' => $stop , 'strand' => $strand);
		}

	# This is a new sequence, so starts a new block by default
	} else {
		
		# finish and print block
                $this_block{'start'} = 1 if ($this_block{'start'} - $max_intron < 0);
                $this_block{'stop'} = ($this_block{'stop'}+$max_intron);
                printf $this_block{'seq_name'} . "\t" . $this_block{'start'} . "\t" . $this_block{'stop'} . "\t" . $this_block{'strand'} . "\n";
		
		%this_block = ( 'seq_name' => $seq, 'start' => $start, 'stop' => $stop, 'strand' => $strand );
	}
	
	$previous_start = $this_start;
	$previous_strand = $strand;
}

# finish and print block
$this_block{'start'} = 1 if ($this_block{'start'} - $max_intron < 0);
$this_block{'stop'} = ($this_block{'stop'}+$max_intron);
printf $this_block{'seq_name'} . "\t" . $this_block{'start'} . "\t" . $this_block{'stop'} . "\t" . $this_block{'strand'} . "\n";
