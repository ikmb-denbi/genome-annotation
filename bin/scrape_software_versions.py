#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'NF-hints': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Blast': ['v_blast.txt', r"blastn: (\S+)"],
    'Hisat2': ['v_hisat2.txt', r"hisat2-align-s version (\S+)"],
    'GenomeThreader': ['v_gth.txt', r"gth \(GenomeThreader\) (\S+)"],
    'RepeatMasker': ['v_rm.txt', r"version (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['NF-hints'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Blast'] = '<span style="color:#999999;\">N/A</span>'
results['Hisat2'] = '<span style="color:#999999;\">N/A</span>'
results['GenomeThreader'] = '<span style="color:#999999;\">N/A</span>'
results['RepeatMasker'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
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
