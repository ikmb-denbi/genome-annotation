#!/usr/bin/env perl
use strict;
use warnings;

=head1 NAME

chunk_chromosome.pl

=head1 SYNOPSIS

chunk_chromosome.pl

=head1 DESCRIPTION

Split a chromosome assembly into chunks (50kb by default).

=head1 OPTIONS

  -fasta_file  the fasta file with chromosome sequences
  -chunk_file  the fasta file of chunks that will be created [default <fasta_file>.chunk]
  -agp_file    the name of the agp file that will be created [default <fasta_file>.agp]
  -size        the size of the chunks, in bases [default 50000]
  -help        displays this documentation with PERLDOC

=cut

use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $help;
my $fasta_file;
my $chunk_file;
my $agp_file;
my $chunk_size;

&GetOptions(
    'fasta_file=s' => \$fasta_file,
    'chunk_file:s' => \$chunk_file,
    'agp_file:s'   => \$agp_file,
    'size:i'       => \$chunk_size,
    'h|help'       => \$help,
) or ($help = 1);

if(! $fasta_file || !-e $fasta_file) {
  print STDERR "Fasta file not given or doesn't exist\n";
}

if ($help) {
  exec('perldoc', $0);
}

$chunk_file = "$fasta_file.chunk" unless $chunk_file;
$agp_file = "$fasta_file.agp" unless $agp_file;
$chunk_size = 50000 unless $chunk_size;

my $seq_in = new Bio::SeqIO(-format=>'Fasta', -file=>$fasta_file);
my $seq_out = new Bio::SeqIO(-format=>'Fasta', -file=>">$chunk_file");
my $assembly;
my $assembly_index = 1;
while ( my $seq = $seq_in->next_seq ) {
    my $name = $seq->id;
    my $chunk_index = 1;
    my $obj_seq = $seq->seq;
    
    while (length($obj_seq) > $chunk_size) {
        my $chunk = substr($obj_seq, 0, $chunk_size, '');
        my $new_seq = new Bio::Seq(-id=>"$name\_$chunk_index", -seq=>$chunk);
        $seq_out->write_seq($new_seq);
        my $obj_start = 1 + (($chunk_index-1) * $chunk_size);
        my $obj_end = $obj_start + $new_seq->length - 1;
        my @assembly = ($name, $obj_start, $obj_end, $assembly_index, 'O', "$name\_$chunk_index", 1, $new_seq->length, '+');
        $assembly .= join("\t", @assembly)."\n";
        $chunk_index++;
        $assembly_index++;
    }
    my $new_seq = new Bio::Seq(-id=>"$name\_$chunk_index", -seq=>$obj_seq);
    $seq_out->write_seq($new_seq);
    my $obj_start = 1 + (($chunk_index-1) * $chunk_size);
    my $obj_end = $obj_start + $new_seq->length - 1;
    my @assembly = ($name, $obj_start, $obj_end, $assembly_index, 'O', "$name\_$chunk_index", 1, $new_seq->length, '+');
    $assembly .= join("\t", @assembly)."\n";
    $assembly_index++;
}

open(AGP, ">$agp_file");
print AGP $assembly;
close(AGP);

