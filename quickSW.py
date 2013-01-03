import os, sys
import numpy as np

match = +2
mismatch = -3
gap = -3
sparse_threshold = -10

score_min_threshold5 = 20
score_min_threshold3 = 20

p5  = 'GGCCGCCTGCAGGAAA'
p5r = 'TTTCCTGCAGGCGGCC'
adap5 = 'ATCTCTCTCAACAACAACAGGCGAAGAGGAAGGAAAGAGAGAGATGGCC'
p3  = 'CGCGCCACCG'
p3r = 'CGGTGGCGCG'
adap3 = 'ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGATCGCG'
Ap3r = 'A'*8 + p3r
p3T  = p3 + 'T'*8

def quickSW(short_seq, long_seq):
    len_s = len(short_seq)
    len_l = len(long_seq)
    M = np.zeros((len_s, len_l))
    M[:] = -9999
    
    compute_count = 0
    
    for j in xrange(len_l): M[0, j] = match if short_seq[0]==long_seq[j] else mismatch
    M[1:, 0] = range(gap, gap*len_s, gap)
    
    for i in xrange(1, len_s):
        for j in xrange(1, len_l):
            M[i, j] = max(M[i-1,j-1]+(match if short_seq[i]==long_seq[j] else mismatch), \
                          M[i,j-1]+gap, \
                          M[i-1,j]+gap)
            compute_count += 1        
    
    #print >> sys.stderr, "compute_count:", compute_count
    return M

def traceback(M, short_seq, long_seq):
    i = M.shape[0] - 1
    j = M.shape[1] - 1
    score = np.max(M[i,:])
    while M[i, j]!=score: j -= 1
    
    mstr1 = '' # for i
    mstr2 = '' # for j
    mstr0 = '' # for |||
    jend = j
    iend = i
    while i > 0 and j > 0:
        if M[i,j] == M[i-1,j-1]+(match if short_seq[i]==long_seq[j] else mismatch):
            mstr1 += short_seq[i]
            mstr2 += long_seq[j]
            mstr0 += '|' if short_seq[i]==long_seq[j] else ' '
            i -= 1
            j -= 1
        elif M[i,j] == M[i,j-1]+gap:
            mstr1 += '-'
            mstr2 += long_seq[j]
            mstr0 += ' '
            j -= 1
        else:
            mstr1 += short_seq[i]
            mstr2 += '-'
            mstr0 += ' '
            i -= 1
            
    print "score:", score
    print mstr2[::-1], j, "-", jend
    print mstr0[::-1]
    print mstr1[::-1], i, "-", iend
    
    return i, iend, j, jend, score

def find_primer(r):
    """
    Either 
    5' ---- polyA + 3'
    polyT + 3'R ---- 5'R
    """
    pass

def find_5prime(r):
    """
    Find 5' either from the beginning or 5'R at the end
    ONLY look at the first 100 bp
    """
    pass
    
    
    

import pdb
def find_missed_adapters(r, f, logf):
    """
    r --- a Bio.SeqRecord
    ASSUME input is from sequencing 2-3k library or > 3k library 
    (hence only finding adapter once at most in the middle)
    
    Also ONLY look at seqIDs that starts at <movie>/<hole>/<start>_<end>
    where start < 500
    """
    offset = 500
    offset_tail = 200
    seq = r.seq.tostring()[offset:-offset_tail] # only find ones that are in the middle
    if len(seq) == 0:
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        return False
    
    headid = r.id[:r.id.rfind('/')]
    s = int(r.id[r.id.rfind('/')+1:r.id.rfind('_')])
    if s > 500:
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        return False
    
    print >> sys.stderr, r.id
    
    """
    look for polyA + 3 --- adapter --- 3'R + polyT
    
    if matches, output segment as: 0 - (jend1+offset), (jend1+offset+j2) - <end> 
    """
    polyA_i = seq.find('A'*6)
    if polyA_i > 0 and seq.find('T'*6, polyA_i, polyA_i+300) > 0: 
        oldseq = seq
        oldoffset = offset   
        seq = seq[polyA_i:]
        offset += polyA_i 
        M = quickSW(Ap3r, seq)
        if np.max(M[len(Ap3r)-1,:]) > score_min_threshold3:
            i1, iend1, j1, jend1, score1 = traceback(M, Ap3r, seq)
            M = quickSW(p3T, seq[jend1:jend1+300])
            #pdb.set_trace()
            if np.max(M[len(p3T)-1,:]) > score_min_threshold3:
                i2, iend2, j2, jend2, score2 = traceback(M, p3T, seq[jend1:jend1+300])
                logf.write("goal3'! {4} split into {0}-{1}, {2}-{3}\n".format(0, jend1+offset, jend1+offset+j2, len(r.seq), headid))   
                #pdb.set_trace()
                f.write(">{0}/0_{1}\n{2}\n".format(headid, jend1+offset, r.seq[:jend1+offset].tostring()))
                f.write(">{0}/{1}_{2}\n{3}\n".format(headid, jend1+offset+j2, len(r.seq), r.seq[jend1+offset+j2:])) 
                return True
        seq = oldseq
        offset = oldoffset   
         
    """
    look for 5'R --- adapter ---- 5'
    
    if matches, output segment as: 0 - (jend1+offset), (jend1+offset+j2) - <end> 
    """
    p5r_i = max(seq.find('TTCCTG'), seq.find('GCAGGC'), seq.find('GCGGCC'), seq.find('CTGCAG'))
    if p5r_i > 0 and max(seq.find('GGAAA',p5r_i,p5r_i+300), seq.find('GCCTGC',p5r_i,p5r_i+300), seq.find('TGCAGG',p5r_i,p5r_i+300), seq.find('GGCCGC',p5r_i,p5r_i+300)) > 0:
        seq = seq[p5r_i-20:]
        offset += p5r_i-20
        M = quickSW(p5r, seq)
        if np.max(M[len(p5r)-1,:]) > score_min_threshold5:
            i1, iend1, j1, jend1, score1 = traceback(M, p5r, seq)
            M = quickSW(p5, seq[jend1:jend1+300])
            #pdb.set_trace()
            if np.max(M[len(p5)-1,:]) > score_min_threshold5:
                i2, iend2, j2, jend2, score2 = traceback(M, p5, seq[jend1:jend1+300])
                logf.write("goal5'! {4} split into {0}-{1}, {2}-{3}\n".format(0, jend1+offset, jend1+offset+j2, len(r.seq), headid))  
                #pdb.set_trace()
                f.write(">{0}/0_{1}\n{2}\n".format(headid, jend1+offset, r.seq[:jend1+offset].tostring()))
                f.write(">{0}/{1}_{2}\n{3}\n".format(headid, jend1+offset+j2, len(r.seq), r.seq[jend1+offset+j2:]))
                return True
     
    # for writing things that did not get split       
    f.write(">{0}\n{1}\n".format(r.id, r.seq))
    return False

if __name__ == "__main__":
    from Bio import SeqIO
    
    f = open('filtered_subreads.asym_adap_find.fasta', 'w')
    logf = open('filtered_subreads.asym_adap_find.log', 'w')
    
    with open('filtered_subreads.fasta') as h:
        for r in SeqIO.parse(h, 'fasta'):
            find_missed_adapters(r, f, logf)
            
    f.close()
    logf.close()
        