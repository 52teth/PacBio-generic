#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

inputs = sys.argv[1].split(",")
output = sys.argv[2] + '.seqlengths.png'
range_min = int(sys.argv[3])
range_max = int(sys.argv[4])

raw = []
for input in inputs:
    raw += [len(r.seq) for r in SeqIO.parse(open(input), 'fasta')]
raw = filter(lambda x: range_min<=x<=range_max, raw)

seqlengths = np.array(raw)

bins = (max(seqlengths)-min(seqlengths))/100 + 1
y,binEdges = np.histogram(seqlengths, bins=bins)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

xnew = np.linspace(bincenters.min(), bincenters.max(), 300)
ysmooth = spline(bincenters, y, xnew)
plt.plot(xnew, ysmooth, '-')
plt.xlabel("Sequence Length")
plt.ylabel("Count")
plt.show()
plt.savefig(output)
