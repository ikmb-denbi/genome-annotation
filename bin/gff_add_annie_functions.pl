#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The GFF name of the file to read. 
    [--annie filename]
		The annie annotation report to read.
  Ouput:    
        Prints to standard output
};

my $annie = undef;
my $gff = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "annie=s" => \$annie);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

open (my $GFF, '<', $gff) or die "FATAL: Can't open file: $gff for reading.\n$!\n";
open (my $ANNIE, '<', $annie) or die "FATAL: Can't open file: $annie for reading.\n$!\n";

# Hold functional annotation data, grouped by transcript model id and functional feature 
my %annie_data ; 

# PARSE ANNIE REPORT
while (<$ANNIE>) {
	
	chomp;
	my $line = $_;
	my ($acc,$feat,$val) = split("\t", $line);
	
	unless (exists $annie_data{$acc}) {
		my %data ;
	   	$annie_data{$acc} = \%data ; 
	}
	
	if ($feat eq "Dbxref") {
		if (exists $annie_data{$acc}) {
			push ( @{ $annie_data{$acc}{$feat} }, $val ) ;
		} else {
			my @data = ( $val ) ;
			$annie_data{$acc}{$feat} = \@data ; 
		}
	}  else {
		$annie_data{$acc}{$feat} = $val ; 
	}	
}

close $ANNIE;

# PARSE GFF FILE
while (<$GFF>) {
	chomp; 
	my $line = $_; 

	my %entry = parse_gff($line);
	
	my $feature_id = $entry{"attributes"}{"ID"} ;
	
	if ( exists $annie_data{$feature_id} ) {
		
		my $annotations = $annie_data{$entry{"attributes"}{"ID"}} ;
		
		if (defined $annotations) {
			
			if ( $entry{"feature"} eq "gene" ) {
				$entry{"attributes"}{"Name"} =  $annotations->{"name"};
			} elsif ( $entry{"feature"} eq "mRNA" || $entry{"feature"} eq "transcript") {
				if ($annotations->{"Dbxref"}) {
					my @go_terms;
					foreach my $xref (@{$annotations->{"Dbxref"}}) {
						if ($xref =~ /GO\:.*/) {
							push(@go_terms, $xref);
						}
					}
					$entry{"attributes"}{"Ontology_term"} = join(",", @go_terms );
				}
				if ($annotations->{"product"}) {
					$entry{"attributes"}{"description"} = $annotations->{"product"};
				}
			}		
		} # defined annotations
		
		gff_print(%entry) ;

	} # exists annie data
}

close $GFF;

# -----------
# - METHODS - 

sub gff_print() {
	
	chomp;
	my %data = @_;
	
	my $attribs = "";
	foreach my $key (keys %{$data{"attributes"}} ) {
		my $value = $data{"attributes"}{$key};
		$attribs .= $key . "=" . $value . ";" ;
	}
	printf $data{"seq_name"} . "\t" . $data{"source"} . "\t" . $data{"feature"} . "\t" . $data{"seq_start"} . "\t" . $data{"seq_end"} . "\t" . $data{"score"} . "\t" . $data{"strand"} . "\t" . $data{"phase"} . "\t" . $attribs ."\n" ;
		
}

sub parse_gff() {
	
	chomp;
	my %data = $_;
	
	my $line = $_;
	my %answer;
	my %attributes;
	
	my @temp = split("\t", $line);
	
	$answer{"seq_name"} = @temp[0];
	$answer{"source"} = @temp[1];
	$answer{"feature"} = @temp[2];
	$answer{"seq_start"} = @temp[3];
	$answer{"seq_end"} = @temp[4];
	$answer{"score"} = @temp[5];
	$answer{"strand"} = @temp[6];
	$answer{"phase"} = @temp[7];
	
	my @attribs = split(";", @temp[8]);
	foreach my $attrib (@attribs) {
		$attrib =~ s/^\s+//;
		my ($key,$value) = split("=", $attrib);
		$attributes{$key} = $value;
	}
	
	$answer{"attributes"} = \%attributes;
	
	return %answer ;
}



