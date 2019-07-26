#!/usr/bin/env python
import sys
from src.annotation import Annotation

#this functions takes an ipr file and returns a list of annotations. 2 types of annotations are retrieved based on the following keys: "name" and "product"
def read_sprot(blast_file, gff_file, fasta_file):
    #retrieve relevant information from files
    fasta_info = get_fasta_info(fasta_file)
    gff_info = get_gff_info(gff_file)
    blast_info = get_blast_info(blast_file)

    sprot_list = []
    for mrna, dbxref in blast_info.items(): #blast_info maps mrna's to dbxrefs
        if dbxref not in fasta_info: #these two if's shouldn't occur but just in case...
            print(mrna+" has dbxref "+dbxref+" that's not in the fasta. Skipping...")
            continue
        if mrna not in gff_info:
            print( mrna+" not in gff. Skipping...")
            continue

        #fasta_info maps dbxrefs to products and names
        product = fasta_info[dbxref][0] 
        gene_name = fasta_info[dbxref][1]
        #gff_info maps mrna's to the parent gene id's
        gene_id = gff_info[mrna]

        #add annotations to annotation list
        sprot_list.append(Annotation(gene_id, "name", gene_name))
        sprot_list.append(Annotation(mrna, "product", product))
    return sprot_list

#this function reads a fasta file and returns a dictionary mapping dbxrefs to the 2-tuple (product, gene name)
def get_fasta_info(fasta_file):
    dbxrefs = {}
    for line in fasta_file:
        if line[0] == '>': #if we have a header line

            words = line.split(" ") #we break the line by spaces to get "words"

            ref = words[0][1:] #we are assuming the first "word" is the dbxref with no spaces. the "[1:]" is to get rid of the starting '>'

            
            i=0
            #loop through the words till we find "OS=" (which we assume exists). We are assuming all the words between here and the ref is the product.
            while words[i].find("OS=") == -1: 
                i += 1
            product = " ".join(words[1:i])
            #loop through the words till we find "GN=" or "PE=". We are assuming "PE=" comes immediately after "GN=" so if we hit "PE=" first, then the gene name doesn't exist. We also assume the gene name is one word
            while (words[i].find("GN=") == -1 and words[i].find("PE=") == -1) and (i+1 < len(words)):
                i += 1
            if not words[i].find("GN=") == -1: #if gene name exists
                name = words[i][3:] #the "[3:]" is to get rid of the "GN=" in the beginning
            else: #if gene name doesn't exist, use part of the dbxref for the name
                name = ref.split("|")[2].split("_")[0]
            #add to dictionary
            dbxrefs[ref] = (product,name)
    return dbxrefs

#this function reads a blast file and returns a dictionary mapping mrna's to dbxrefs. This function assumes the file is tab-separated and has mrna in the 0th and dbxref in the 1st column
def get_blast_info(blast_file):
    mrna_dbxrefs = {}
    for line in blast_file:
        columns = line.split("\t")
        mrna = columns[0]
        ref = columns[1]
        if mrna not in mrna_dbxrefs:
            mrna_dbxrefs[mrna] = ref
    return mrna_dbxrefs

#this function reads a gff file and returns a dictionary mapping mrna id's to its parent gene id
def get_gff_info(gff_file):
    mrna_genes = {}
    for i, line in enumerate(gff_file):
        columns = line.split("\t")
        if len(columns)>1 and columns[2] == "mRNA": #if this is an "mRNA row"
            mrna_id = ""
            parent_gene = ""
            for attribute in columns[8].strip().split(";"):
                split = attribute.split("=")
                key, val = split[0], split[1]
                if key == "ID":
                    mrna_id = val
                elif key == "Parent":
                    parent_gene = val
            if not mrna_id or not parent_gene:
                print("Failed to get mRNA info at line "+str(i)+" of GFF because it is missing the ID and/or Parent attributes")
                continue
            mrna_genes[mrna_id] = parent_gene
    return mrna_genes


