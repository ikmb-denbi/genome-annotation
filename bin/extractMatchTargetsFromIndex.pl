#!/usr/bin/env perl
# Get proteins in FASTA format from a cdb index using a list

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
my $outfile = undef;
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

while (<$IN>) {
	my $line = $_;
	chomp $line;

	my @temp = split(/\t+/,$line);
	my $query = $temp[0];

	my $query_clean = $query;
	$query_clean =~ s/\|/_/g ;
	
	
	my $fa_clean = "$query_clean._query_.fa" ;

	unless (-e $fa_clean) {
		#printf "Will extract $query to $fa_clean\n";
		system("cdbyank $db -a \"${query}\" > $fa_clean");
	}
}
