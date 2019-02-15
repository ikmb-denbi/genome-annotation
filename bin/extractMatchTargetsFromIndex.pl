#!/usr/bin/perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--matches filename]
		The name of the file to read. 
    [--db filename]
		The name of the CDBtools data index

};

my $matches = undef;
my $db = undef;
my $help;

GetOptions(
    "help" => \$help,
    "matches=s" => \$matches,
    "db=s" => \$db,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open (my $IN, '<', $matches) or die "FATAL: Can't open file: $matches for reading.\n$!\n";

while (<$IN) {
	$line = $_;
	chomp $line;

	@temp = split(/\t+/,$line);
	$query = $temp[0];
	$target = $temp[1];

	my $query_clean = $query;
	$query_clean =~ s/\|/_/g ;

	my $fa_clean = $query_clean + ".fa" ;

	unless (-e $fa_clean) {
		system("cdbyank $db -a '$query' > $fa_clean");
	}
}
