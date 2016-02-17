#!/usr/bin/env python
from Bio import SeqIO
import sys
file = sys.argv[1]
ids = sys.argv[2:]

for r in SeqIO.parse(open(file), 'fasta'):
    if any(r.id.startswith(x) for x in ids):
        print ">"+r.id
        print r.seq

