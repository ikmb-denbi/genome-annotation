#!/bin/env ruby

require 'bio'
require 'optparse'
require 'ostruct'

module GFF

	class Entry
		attr_accessor :seq, :source, :feature, :start, :stop, :score, :strand, :phase, :attributes

		def initialize(line)
			@seq,@source,@feature,@start,@stop,@score,@strand,@phase,tmp = line.split("\t").collect{|f| f.strip }
			@start = @start.to_i
			@stop = @stop.to_i
			@attributes = _make_attributes(tmp)	
		end

		def _make_attributes(line)
		  puts self.inspect if line.nil?
            
			answer = {}
			elements = line.strip.split(";")
			elements.each { |e| answer[e.split("=")[0]] = e.split("=")[1] }
			return answer
		end

		def output
			attributes = @attributes.collect{|k,v| "#{k}=#{v}" }.join(";")
			return "#{@seq}\t#{@source}\t#{@feature}\t#{@start}\t#{@stop}\t#{@score}\t#{@strand}\t#{@phase}\t#{attributes}"
		end
		
		def strip_attributes!
		  self.attributes.delete("inference")
      self.attributes.delete("product")
      self.attributes.delete("locus_tag")
      self.attributes.delete("gene")
		end
		
	end

end

def parse_annie(lines)
  
  bucket = {}
  
  lines.each do |line|
    
    accession,feature,value = line.strip.split("\t").collect{|e| e.strip}
        
    unless bucket.has_key?(accession)
      bucket[accession] = {}
    end
    
    if feature == "Dbxref"
      bucket[accession].has_key?(feature) ? bucket[accession][feature] << value : bucket[accession][feature] = [ value ]
    else
      bucket[accession][feature] = value
    end
    
  end
  
  return bucket
  
end

options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-g","--gff", "=GFF","GFF annotation") {|argument| options.gff = argument }
opts.on("-a","--annie", "=ANNIE","Annie annotation") {|argument| options.annie = argument }
opts.on("-h","--help","Display the usage information") {
  puts opts
  exit
}

opts.parse!

annie = parse_annie(IO.readlines(options.annie))


gff = File.open(options.gff,"r")

while (line = gff.gets)
    
  if line.match(/^##.*/)
      puts line
  elsif line.match(/^#.*/)
    # do nothing, i.e. remove these lines. 
  else
        
    entry = GFF::Entry.new(line)
    
    if entry.feature == "gene"
     
     if annie.has_key?(entry.attributes["ID"])
       xref = annie[entry.attributes["ID"]]
       entry.attributes["Name"] = xref["name"]
     end
     
    elsif entry.feature == "mRNA" or entry.feature == "transcript"
      if annie.has_key?(entry.attributes["ID"])
        xrefs = annie[entry.attributes["ID"]] 
        if xrefs.has_key?("Dbxref")
          entry.attributes["Dbxref"] = xrefs["Dbxref"].select{|x| x.include?("GO:") == false }.join(",").gsub(/\|/, ',')
          go_terms = xrefs["Dbxref"].select{|x| x.include?("GO:") }
          unless go_terms.empty?
            entry.attributes['Ontology_term'] = go_terms.join(",").gsub(/\|/, ',')
          end
        end
        if xrefs.has_key?("product")
          entry.attributes["description"] = xrefs["product"]
        end
      end
    end
      
    puts entry.output
  end
    
end

gff.close
