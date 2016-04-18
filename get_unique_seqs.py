#!/usr/bin/env python
import os, sys
from Bio import SeqIO

input = sys.argv[1]
seen = set()

f = open(input[:input.rfind('.')]+'.derep.fasta', 'w')
for r in SeqIO.parse(open(input),'fasta'):
    if r.id in seen: continue
    seen.add(r.id)
    f.write(">{0}\n{1}\n".format(r.id, r.seq))
f.close()
