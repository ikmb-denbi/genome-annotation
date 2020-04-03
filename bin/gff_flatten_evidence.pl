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

my $previous_seq = undef;
my $previous_group = undef;

my %container;

while (<$IN>) {

	my $line = $_;
	chomp $line;

	my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);	

	my %entry = { "seq" => $seq, "source" => $source, "feature" => $feature, "start" => $start, "stop" => $stop, "phase" => $phase, "strand" => $strand, "score" => $score } ;

	my %attribs;
	
	my @fields = split(";",$info);

	foreach my $f (@fields) {
		my ($key,$value) = split("=",$f);
		$attribs{$key} = $value;
	}
	
	$entry{"attributes"} = \%attribs;

	my $this_group = $attribs{"grp"};

	#printf STDERR $this_group . "\n" ;

	if (defined $container{$this_group} ) {
		push(@{ $container{$this_group}{'elements'} },\%entry);
	} else {
		my @data = ( %entry );
		$container{$this_group}{'elements'} = \@data ;
		$container{$this_group}{'seq'} = $seq;
		$container{$this_group}{'start'} = $start;
	}

	$previous_seq = $seq;

}

sub build_clusters {


	my %container = $_ ;
	#printf Dumper(\%container);

	my $count = scalar (keys %container);
        printf STDERR "Would print $count now!\n";
	
}

