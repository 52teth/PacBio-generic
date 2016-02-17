#!/usr/bin/env python
import os, sys
from Bio import SeqIO

fastafile = sys.argv[1]
listfile = sys.argv[2]

seqs = [line.strip() for line in open(listfile)]
for r in SeqIO.parse(open(fastafile), 'fasta'):
    if r.id in seqs: 
        print ">" + r.id
        print r.seq
