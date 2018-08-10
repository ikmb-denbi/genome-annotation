#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os.path

regexes = {
    'NF-hints': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': "not used",
    'Blast': "not used",
    'GenomeThreader': "not used",
    'RepeatMasker': "not used",
    'Trim Galore!': "not used",
    'Hisat2': "not used",
#    'Trinity': "not used",
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}

if os.path.exists("v_fastqc.txt"):
	regexes['FastQC'] = ['v_fastqc.txt', r"FastQC v(\S+)"]

if os.path.exists('v_blast.txt'):
	regexes['Blast'] = ['v_blast.txt', r"blastn: (\S+)"]

if os.path.isfile('v_gth.txt'):
	regexes['GenomeThreader'] = ['v_gth.txt', r"gth \(GenomeThreader\) (\S+)"]

if os.path.isfile('v_rm.txt'):
	regexes['RepeatMasker'] = ['v_rm.txt', r"version (\S+)"]

if os.path.isfile("v_trim_galore.txt"):
	regexes['Trim Galore!'] = ['v_trim_galore.txt', r"version (\S+)"]

if os.path.isfile("v_hisat2.txt"):
	regexes['Hisat2'] = ['v_hisat2.txt', r"hisat2-align-s version (\S+)"]

if os.path.isfile('v_trinity.txt'):
	regexes['Trinity'] = ['v_trinity.txt', r"Trinity version: Trinity-(\S+)"]


results = OrderedDict()
results['NF-hints'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Blast'] = '<span style="color:#999999;\">N/A</span>'
results['GenomeThreader'] = '<span style="color:#999999;\">N/A</span>'
results['RepeatMasker'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['Hisat2'] = '<span style="color:#999999;\">N/A</span>'
#results['Trinity'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
	if v != "not used":
		with open(v[0]) as x:
			versions = x.read()
			match = re.search(v[1], versions)
			if match:
				results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-hints-software-versions'
section_name: 'NF-hints Software Versions'
section_href: 'https://git.ikmb.uni-kiel.de/m.torres/NF-hints.git'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
