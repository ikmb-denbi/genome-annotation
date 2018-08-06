#!/usr/bin/env perl

use strict;
use warnings;

my $repeatfile = $ARGV[0];
my @t;

open(REPEATS, "$repeatfile") || die "ERROR: cannot open $repeatfile\n";

while (<REPEATS>) {
	chomp $_; 
	$_ =~ s/^\s+//; 
	@t = split(/\s+/,$_); 
	print STDOUT $t[4]."\t"."repmask\tnonexonpart\t".$t[5]."\t".$t[6]."\t0\t.\t.\tsrc=RM\n";
}
close REPEATS
