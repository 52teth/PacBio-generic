#!/usr/bin/env python
import os, sys
import GFF
import numpy as np
from Bio import SeqIO

epsilon=20
min_coverage=95

def validate_53seen(gtf, chrom, pos5, pos3, strand):
    """
    Given the 5'-3' position (probably from GMAP)
    See if they fall at the first/last exon of any transcript
    Returns hit5, hit3
    
    For hit5, (i-th exon, distance to start of i-th exon, id) where argmin i
    For hit3, (-j-th exon, distance to end of j-th exon, id) where argmin j
    """
    hit5 = (np.inf, np.inf, 'NA') # i-th exon, dist
    hit3 = (np.inf, np.inf, 'NA') # -i-th exon, dist
    assert pos5 < pos3 # if strand is -, then actually 3'-5'
    for id in gtf.find(chrom, pos5, pos3):
        if gtf.transcript_info[id]['strand']!=strand: continue
        exons = gtf.get_exons(id)
        for i in xrange(len(exons)):
            if exons[i].start-epsilon<=pos5<=exons[i].end+epsilon:
                d5 = abs(pos5-exons[i].start)                
                if hit5[0] > i or (hit5[0] == i and d5 < hit5[1]):
                    hit5 = (i,d5,id)
            j = len(exons)-i-1
            if exons[i].start-epsilon<=pos3<=exons[i].end+epsilon:
                d3 = abs(pos3-exons[i].end)
                if hit3[0] > j or (hit3[0]==j and d3 < hit3[1]):
                    hit3 = (j,d3,id)
    if strand == '+':
        return hit5, hit3
    else:
        return hit3, hit5

def validate_polyA(gtf, chrom, pos5, pos3, strand):
    """
    Given the 3' position
    Find the closest polyA site (+- 100bp)
    """
    dist = np.inf
    if strand == '-':
        pos5, pos3 = pos3, pos5
    for id in gtf.find(chrom, pos3-100, pos3+100):
        if gtf.transcript_info[id]['strand']!=strand: continue
        exons = gtf.get_exons(id)
        assert len(exons)==1
        if strand == '+':
            dist = min(dist, abs(pos3-exons[0].start))
        else:
            dist = min(dist, abs(pos3-exons[0].end))
    return dist


def main(gmap_filename, fasta_filename):
    """
    Given a GMAP output (.gff) compare the aligned start/end
    to Gencode annotations (transcript & polyA)

    Need the original fasta to get sequence length
    """
    seqlen_dict = dict([(r.id,len(r.seq)) for r in SeqIO.parse(open(fasta_filename),'fasta')])
    gtf_f = '/home/UNIXHOME/etseng/share/gencode/gencode.v15.annotation.gtf'
    gtfA_f = '/home/UNIXHOME/etseng/share/gencode/gencode.v15.polyAs.gtf'

    gtf = GFF.GTF(gtf_f)
    gtfA = GFF.polyAGFF(gtfA_f)

    f = open(gmap_filename+'.summary', 'w')
    f.write("ID\thit5_exon\thit5_dist\thit5_id\thit3_exon\thit3_dist\thit3_id\thitA_dist\n")
    reader = GFF.gmapGFFReader(gmap_filename)
    while True:
        try:
            r = reader.next()
        except AssertionError: #ignore bad gmap output
            continue
        except StopIteration:
            break
        except:
            continue
        if r.coverage < min_coverage: continue
        # IMPORTANT! if r.start/r.end is not complete, extend it!
        r_start_corrected = r.start - r.seq_exons[0].start
        r_end_corrected = r.end + (seqlen_dict[r.seqid] - r.seq_exons[-1].end)
        hit5, hit3 = validate_53seen(gtf, r.chr, r.start, r.end, r.strand)
        hitA = validate_polyA(gtfA, r.chr, r.start, r.end, r.strand)
        f.write("{id}\t{e5}\t{d5}\t{i5}\t{e3}\t{d3}\t{i3}\t{dA}\n".format(\
                id=r.seqid, e5=hit5[0], d5=hit5[1], i5=hit5[2], e3=hit3[0], d3=hit3[1], i3=hit3[2], dA=hitA))

    f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
