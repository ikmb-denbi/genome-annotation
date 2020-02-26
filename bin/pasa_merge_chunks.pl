#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Cwd;
use Data::Dumper;

# Take list of primary pasa outputs (assemblies.fasta, pasa_assemblies.gff3) and merge into one file

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:

  Ouput:    

};

my $help;
my $base = "pasa_db_merged" ;

GetOptions(
    "base=s" => \$base,
    "help" => \$help);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $dir = getcwd;

opendir(DIR, $dir) or die $!;

my $file_base = undef;

my $fasta = $base . ".assemblies.fasta";
open(FASTA_OUT,">",$fasta) or die "Couldnt open fasta out\n";

my $gff = $base . ".pasa_assemblies.gff" ;
open(GFF_OUT,">",$gff) or die "Couldnt open gff out\n";

my %bucket_fa;
my %bucket_gff;
my $file_base = undef;

my $align_id = 1;
my $asmbl_id = 1;

# FInd all GFF files
while (my $file = readdir(DIR)) {

        # We only want files
        next unless (-f "$dir/$file");

	my $num = (split /\./,$file)[-4];

	if ($file =~ m/\.pasa_assemblies\.gff3$/) {
			
		$bucket_gff{$num} = $file ;
		
	} elsif ($file =~ m/\.fasta$/) {
		next if ($file =~ $fasta) ;
		$bucket_fa{$num} = $file ;
	}
}

# Find all FASTA files
foreach ( sort {$a<=>$b} keys %bucket_fa) {
        #print "key $_ value: $bucket_fa{$_}\n";

	my $num = $_;
	my $file = $bucket_fa{$num};

	my $fh = IO::File->new();
	$fh->open( $file);
	
	while (<$fh>) {
	
		chomp;
		my $line = $_ ;
		
	
		if ($line =~ /^>.*/ ) {
			print FASTA_OUT ">asmbl_$asmbl_id\n";
			$asmbl_id += 1;
		} else {
			print FASTA_OUT $line . "\n";
		}
	}
	close($fh);
}

close(FASTA_OUT);

# Running ID for assemblies
my $asmbl_id = 1;

# Holding context-specific asmbl ids to detect switches to a new id between loci
my $original_asmbl_id = undef;
my $original_align_id = undef;

# Iterate over all GFF files, sorted by Nextflow counter so they are in the right order
foreach ( sort {$a<=>$b} keys %bucket_gff) {

	my $num = $_;
	my $file = $bucket_gff{$num};
	
	my $fh = IO::File->new();
        $fh->open( $file);
	
	while (<$fh>) {

                chomp;
                my $line = $_ ;

		if ($line =~ /^#.*/) {
			next;
		}

		my ($seq,$source,$feature,$start,$stop,$phase,$strand,$score,$info) = split("\t",$line);

		$source = "PASA_transcript_assemblies" ;

                my @attributes ;

		# Split the GFF attributes
                foreach my $attrib (split ";",$info) {

                        my ($key,$value) = split("=",$attrib);

			# the cDNA alignment id
                        if ($key eq "ID" || $key eq "Parent") {

				if (!defined $original_align_id) {
                                        $original_align_id = $value;
					$value = "align_$align_id";
                                } elsif ( $original_align_id eq $value) {
                                        $value = "align_$align_id";
                                } else {
					$original_align_id = $value;
					$align_id += 1;
					$value = "align_$align_id";
                                }
			
			# The assigned pasa locus/assembly id
                        } elsif ($key eq "Target") {
				
				my @values = split(" ", $value);
				my $tmp  = shift @values ;
				my $id = undef;

	                        printf "original: $original_asmbl_id This: $tmp\n";
				if (!defined $original_asmbl_id) {
					#printf "No id defined, starting new one\n";
                                        $original_asmbl_id = $tmp;
                                        $id = "asmbl_$asmbl_id";
                                } elsif ( $original_asmbl_id eq $tmp) {
					#printf "Is the same asmbl id, not changing\n";
                                        $id = "asmbl_$asmbl_id";
                                } else {
					#printf "Is a different asmbl id, upping counter\n";
                                        $original_asmbl_id = $tmp;
                                        $id = "asmbl_$asmbl_id";
                                        $asmbl_id += 1;
                                }

				my @elements = split(" ", $value);
				my $id_ref = shift @elements;
				my $new_id_ref = "asmbl_" . $asmbl_id;
				my $trunk = join " ",@elements ;
				$value = $new_id_ref . " " . $trunk ;
				
			}
                        push(@attributes,"$key=$value");
                }

		my $result = $seq . "\t" . $source . "\t" . $feature . "\t" . $start . "\t" . $stop . "\t" . $phase . "\t" . $strand . "\t" . $score . "\t" . join(";",@attributes) . "\n";		
		print GFF_OUT $result ;
        }
        close($fh);
}

close(GFF_OUT);

closedir(DIR);
exit 0;
