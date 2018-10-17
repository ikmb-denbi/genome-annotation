#!/bin/env ruby
# == NAME
# script_skeleton.rb
#
# == USAGE
# ./augustus_add_exons.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to fix AUGUSTUS gff files
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-i","--infile", "=INFILE","Input file") {|argument| options.infile = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

options.infile ? input_stream = File.open(options.infile,"r") : input_stream = $stdin
options.outfile ? output_stream = File.new(options.outfile,'w') : output_stream = $stdout

puts "##gff-version 3"

gene_id = nil
transcript_id = nil
cds_counter = 0

while (line = input_stream.gets)
  
#  next if line.match(/^#.*/)
  
  elements = line.strip.split("\t")
  
  feature = line.strip.split("\t")[2]
  
  if feature == "gene" 
    
    gene_id = elements[-1]
    
    elements[-1] = "ID=#{gene_id}"
    
    output_stream.puts elements.join("\t")
    
  elsif feature == "mRNA"
    
    cds_counter = 0
    
    elements[2] = "mRNA"
    
    transcript_id = elements[-1]
    
    elements[-1] = "ID=#{transcript_id};Parent=#{gene_id}"   
    
    output_stream.puts elements.join("\t")
    
  elsif feature == "CDS"
    
    cds_counter += 1
    
    elements[-1] = "ID=#{transcript_id}.CDS-#{cds_counter};Parent=#{transcript_id}"
    
    cds_line = elements.join("\t")
    
    output_stream.puts cds_line
    
    elements[-1] = "ID=#{transcript_id}.EXON-#{cds_counter};Parent=#{transcript_id}"
    
    elements[2] = "exon"
    elements[7] = "."
    
    output_stream.puts elements.join("\t")
    
  end
  
end

input_stream.close
output_stream.close
