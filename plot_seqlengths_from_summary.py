#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from csv import DictReader

input = sys.argv[1]
output = sys.argv[2] + '.seqlengths.png'

raw = []
#Movie,ReadId,#Bases,Readlength,ReadScore,SequencingZMW,Productivity,PassedFilter
for r in DictReader(open(input)):
    if r['PassedFilter']=='1':
        assert r['SequencingZMW']=='1' and r['Productivity']=='1'
        raw.append(int(r['Readlength']))
raw = filter(lambda x: 0<x<15001, raw)

seqlengths = np.array(raw)


bins = (max(seqlengths)-min(seqlengths))/400 + 1
plt.hist(seqlengths, bins=bins, color="#9966FF")
plt.xlabel("Sequence Length")
plt.ylabel("Count")
plt.show()
plt.savefig(output)
