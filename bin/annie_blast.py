#!/usr/bin/env python
import os
import sys
import argparse
from src.sprot import read_sprot
from src.sprot import get_fasta_info
from src.sprot import get_blast_info
from src.sprot import get_gff_info
from src.annotation import write_annotations
#from src.fix import fix_anno

def main(args):
    parser = argparse.ArgumentParser(
    epilog="""
    Modified from annie.py
    Docs at http://genomeannotation.github.io/annie/
    Bugs and feature requests at https://github.com/genomeannotation/annie/issues
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-b', '--blast-output')
    parser.add_argument('-g', '--gff', help="GFF3 file corresponding to assembly")
    parser.add_argument('-db', '--blast-database', help="The fasta file against which BLAST was run")
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    # Make sure we got enough args
    sprot = False
    if args.blast_output and args.gff and args.blast_database:
        sprot = True
    if not (sprot):
        sys.stderr.write("Error: must provide --blast-output, --gff and --blast-database\n\n")
        parser.print_help()
        sys.exit()

    # Open output file
    out = "annie_output.tsv"
    if args.output:
        out = args.output
    outfile = open(out, 'w')

    # Create an empty list to store Annotation objects
    annotations = []

    # Add SwissProt results if requested
    if sprot:
        try:
            blast_file = open(args.blast_output, 'r')
            gff_file = open(args.gff, 'r')
            fasta_file = open(args.blast_database, 'r')
        except IOError:
            print("Sorry, Annie says either one of the files doesn't exist or it could not be read.")
            exit()
        annotations.extend(read_sprot(blast_file, gff_file, fasta_file))
        blast_file.close()
        gff_file.close()
        fasta_file.close()

    #write the annotations to file and close
    write_annotations(annotations, outfile)
    outfile.close()

####################################################################################################

if __name__ == "__main__":
    main(sys.argv)
