#!/bin/env perl
use strict;
use warnings;

my %PATH_COUNTER;

my $bam = $ARGV[0] ;

open BAM,"samtools view $bam |";

        while(<BAM>){
        next if(/^(\@)/);  ## skipping the header lines (if you used -h in the samools command)
        s/\n//;  s/\r//;  ## removing new line
        my @sam = split(/\t+/);  ## splitting SAM line into array

	my %entry = ( "qname" => $sam[0], "flag" => $sam[1], "rname" => $sam[2], "pos" => $sam[3], "mapq" => $sam[4], 
		"cigar" => $sam[5], "rnext" => $sam[6], "pnext" => $sam[7], "tlen" => $sam[8], "seq" => $sam[9], "qual" => $sam[10] );

	my $num_mismatches = 0;

	if ( join("\t",@sam) =~ /NM:i:(\d+)/) {
		$num_mismatches = $1;
        }

	my $strand = ($entry{"flag"} == 0) ? "+" : "-" ;
	my $read_name = $entry{"qname"} ;
	my $scaff_name = $entry{"rname"};
	my ($genome_coords_aref, $query_coords_aref) = get_aligned_coords(%entry);

	my $align_len = 0;

	foreach my $coordset (@$genome_coords_aref) {
                $align_len += abs($coordset->[1] - $coordset->[0]) + 1;
        }

	my $per_id = sprintf("%.1f", 100 - $num_mismatches/$align_len * 100); 

	my $align_counter = "$read_name.p" . ++$PATH_COUNTER{$read_name};

	my @genome_n_trans_coords;
        
        while (@$genome_coords_aref) {
            my $genome_coordset_aref = shift @$genome_coords_aref;
            my $trans_coordset_aref = shift @$query_coords_aref;

            my ($genome_lend, $genome_rend) = @$genome_coordset_aref;
            my ($trans_lend, $trans_rend) = sort {$a<=>$b} @$trans_coordset_aref;

            push (@genome_n_trans_coords, [ $genome_lend, $genome_rend, $trans_lend, $trans_rend ] );

        }

	my @merged_coords;
        push (@merged_coords, shift @genome_n_trans_coords);

        my $MERGE_DIST = 10;
        while (@genome_n_trans_coords) {
            my $coordset_ref = shift @genome_n_trans_coords;
            my $last_coordset_ref = $merged_coords[$#merged_coords];
            
            if ($coordset_ref->[0] - $last_coordset_ref->[1] <= $MERGE_DIST) {
                # merge it.
                $last_coordset_ref->[1] = $coordset_ref->[1];

                if ($strand eq "+") {
                    $last_coordset_ref->[3] = $coordset_ref->[3];
                } else {
                    $last_coordset_ref->[2] = $coordset_ref->[2];
                }
            }
            else {
                # not merging.
                push (@merged_coords, $coordset_ref);
            }
        }

	foreach my $coordset_ref (@merged_coords) {
            my ($genome_lend, $genome_rend, $trans_lend, $trans_rend) = @$coordset_ref;
            
            print join("\t",
                       $scaff_name,
                       "genome",
                       "cDNA_match",
                       $genome_lend, $genome_rend,
                       $per_id,
                       $strand,
                       ".",
                       "ID=$align_counter;Target=$read_name $trans_lend $trans_rend") . "\n";
        }
        print "\n";
        
        
}


sub get_aligned_coords {

	my %entry = @_;

	my $genome_lend = $entry{"pos"};
	my $alignment = $entry{"cigar"};
	my $query_lend = 0;

	my @genome_coords;
	my @query_coords;

	$genome_lend--;

	while ($alignment =~ /(\d+)([A-Z])/g) {

		my $len = $1;
		my $code = $2;
		
		unless ($code =~ /^[MSDNIH]$/) {
			exit 1;  "Error, cannot parse cigar code [$code] ";
		}
		
		# print "parsed $len,$code\n";
		
		if ($code eq 'M') { # aligned bases match or mismatch
			
			my $genome_rend = $genome_lend + $len;
			my $query_rend = $query_lend + $len;
			
			push (@genome_coords, [$genome_lend+1, $genome_rend]);
			push (@query_coords, [$query_lend+1, $query_rend]);
			
			# reset coord pointers
			$genome_lend = $genome_rend;
			$query_lend = $query_rend;
			
		}
		elsif ($code eq 'D' || $code eq 'N') { # insertion in the genome or gap in query (intron, perhaps)
			$genome_lend += $len;
			
		}

		elsif ($code eq 'I'  # gap in genome or insertion in query 
               ||
               $code eq 'S' || $code eq 'H')  # masked region of query
        { 
            $query_lend += $len;

		}
	}

	return(\@genome_coords, \@query_coords);
}
