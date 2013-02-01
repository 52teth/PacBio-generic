#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

input = sys.argv[1]
output = sys.argv[1] + '.seqlengths.png'

raw = [len(r.seq) for r in SeqIO.parse(open(input), 'fasta')]
#raw = filter(lambda x: 0<x<3001, raw)

seqlengths = np.array(raw)


bins = (max(seqlengths)-min(seqlengths))/100 + 1
plt.hist(seqlengths, bins=bins)
plt.xlabel("Sequence Length")
plt.ylabel("Count")
plt.show()
plt.savefig(output)
