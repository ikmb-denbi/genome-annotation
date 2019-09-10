#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--annie filename]
		The annie annotation report to read.
    [--gff filename]
                The GFF name of the file to read.
    [--fasta filenames]
		The FASTA name of the file to read.
  Ouput:    
    [--out_gff]
		The GFF output file.
    [--out_fasta]
		The FASTA output file.
};

my $annie = undef;
my $gff = undef;
my $fasta = undef;
my $out_gff = undef;
my $out_fasta = undef;
my $help;

GetOptions(
    "help" => \$help,
    "annie=s" => \$annie,
    "gff=s" => \$gff,
    "fasta=s" => \$fasta,
    "out_gff=s" => \$out_gff,
    "out_fasta=s" => \$out_fasta);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

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

open (my $GFF, '<', $gff) or die "FATAL: Can't open file: $gff for reading.\n$!\n";
open (my $OUTGFF, ">$out_gff") or die "FATAL: Can't open file: $out_gff for writing.\n";

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
close $OUTGFF;

open (my $FASTA, '<', $fasta) or die "FATAL: Can't open file: $fasta for reading.\n$!\n";
open (my $OUTFASTA, ">$out_fasta") or die "FATAL: Can't open file: $out_fasta for writing.\n";

while (<$FASTA>) {
        chomp;
	my $line = $_;
	
	if ($line =~ /^>/) {
		my @words = split(" ",$line);
		my $gene_ID = $words[0]; # Get only the gene ID (first word)
		$gene_ID =~ s/^.//s; # Remove the ">" sign 
		
		if ( exists $annie_data{$gene_ID} ) {

			my $annotation = $annie_data{$gene_ID};

			if ($annotation->{"product"}) { # In most cases, the first word in the header will be the transcript ID
				my $description = $annotation->{"product"};

			} elsif ($annotation->{"name"}) { # It could be that the first word in the header is the gene ID. I am also interested in the gene ID
				my $description = $annotation->{"name"};
			}
			# Description will be gene ID if the header has gene ID first, or gene Name if the header has transcript ID first
			printf $OUTFASTA $line." Description=".$description."\n"; 

		} else {

			printf $OUTFASTA $line." Description=Unknown\n";
		}
	} else { # If it's not a header but just sequence:
		printf $OUTFASTA "$line\n";
	}
}

close $FASTA;
close $OUTFASTA;

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
	printf $OUTGFF $data{"seq_name"} . "\t" . $data{"source"} . "\t" . $data{"feature"} . "\t" . $data{"seq_start"} . "\t" . $data{"seq_end"} . "\t" . $data{"score"} . "\t" . $data{"strand"} . "\t" . $data{"phase"} . "\t" . $attribs ."\n" ;
		
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



