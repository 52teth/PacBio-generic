#!/usr/bin/env  python
import os, sys, numpy as np
from Bio import SeqIO

file = sys.argv[1]
lens = np.array([len(r.seq) for r in SeqIO.parse(open(file), 'fasta')])

print "Number of seqs:", len(lens)
print "Min:", np.min(lens)
print "Max:", np.max(lens)
print "Median:", np.median(lens)
print "Mean:", np.mean(lens)

