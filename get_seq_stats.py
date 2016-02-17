#!/usr/bin/env python
import os, sys
import numpy as np
from Bio import SeqIO

file = sys.argv[1]

f = open(file + '.seqlengths.txt', 'w')
lens = []
for r in SeqIO.parse(open(file),'fasta'):
    f.write(r.id + '\t' + str(len(r.seq)) + '\n')
    lens.append(len(r.seq))
f.close()

lens = np.array(lens)

print "{0} sequences".format(len(lens))
print "min:", np.min(lens)
print "max:", np.max(lens)
print "avg:", np.mean(lens)

# print by 1 kb bins
for i in xrange(0, np.max(lens)/1000+1):
    print "{0}-{1} kb: {2}".format(i, i+1, sum((i*1000 <= lens) & (lens < (i+1)*1000)))

