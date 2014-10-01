#!/usr/bin/env python
import os
from Bio import SeqIO
from csv import DictReader
seqlen_dict=dict((r.id,len(r.seq)) for r in SeqIO.parse(open('reads_of_insert.fasta'),'fasta'))
bins=[0]*24 # every 500 bp
from collections import defaultdict
bins=defaultdict(lambda: {'fl':0, 'nofl':0})
reader=DictReader(open('isoseq_primer_info.csv'),delimiter=',')
# header: id      strand  fiveseen        polyAseen       threeseen       fiveend polyAend        threeend        primer       chimera
for r in reader:
    b=seqlen_dict[r['id'][:r['id'].rfind('/')]+'/ccs']/500
    if r['fiveseen']=='1' and r['threeseen']=='1' and r['polyAseen']=='1': bins[b]['fl']+=1
    else: bins[b]['nofl']+=1
    
f=open('isoseq_flnc.fasta.bin_by_len.txt','w')
dirname=os.path.abspath('.').split('/')[-1]
f.write(dirname+',')
for i in xrange(24): f.write(str(bins[i]['fl']*100./max(1, bins[i]['fl']+bins[i]['nofl']))+',')
f.write('\n')
