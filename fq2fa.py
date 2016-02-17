#!/usr/bin/env python
from Bio import SeqIO
import sys
file = sys.argv[1]
if file.endswith('.fq'):
    out = file[:-2] + 'fa'
elif file.endswith('.fastq'):
    out = file[:-5] + 'fasta'
else:
    raise Exception, "{0} not .fa or .fasta".format(file)

with open(out, 'w') as f:
    for r in SeqIO.parse(open(file),'fastq'):
        f.write(">{0}\n{1}\n".format(r.id, r.seq))

