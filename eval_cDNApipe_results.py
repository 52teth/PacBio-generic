import os, sys
import numpy as np
from collections import defaultdict
from Bio import SeqIO

import miscBio
import eval_accuracy_by_gmap
import eval_refmap_by_blasr

def get_zmw(filename):
    """
    Returns dict of {zmw} --> count of sequences (for CCS, must be exactly 1)
    """
    d = defaultdict(lambda: 0)
    for r in SeqIO.parse(open(filename), 'fasta'):
        if r.id.endswith('/ccs'): #  old version of CCS id
            zmw = r.id[:-4]
        elif r.id.count('/') == 1: # new version of CCS
            zmw = r.id
        elif r.id.count('/') == 2: # subread
            zmw = r.id[:r.id.rfind('/')]
        else:
            raise Exception, "Unrecognized ID format: {0}".format(r.id)
        d[zmw] += 1
    return d

def get_readlen(filename): return dict((r.id, len(r.seq)) for r in SeqIO.parse(open(filename), 'fasta'))

def summarize_CCS(ccs_filename='ccs_reads.fasta', nonccs_filename='nonccs_subreads.fasta'):
    """
    1. % of CCS (per total number of ZMW)
    2. Number of CCS
    3. Number of non-CCS subread
    4. Mean CCS readlength
    5. Mean non-CCS subread readlength
    """
    d_ccs = get_zmw(ccs_filename)
    d_nonccs = get_zmw(nonccs_filename)
    
    r_ccs = get_readlen(ccs_filename)
    r_nonccs = get_readlen(nonccs_filename)
    
    a = len(d_ccs)
    b = len(d_nonccs) + len(d_ccs)
    print("% of ZMWs with CCS: {0}/{1} ({2:.1f}%)".format(a, b, a*100./b))
    print("Number of CCS reads: {0}".format(len(d_ccs)))
    print("Number of non-CCS subreads: {0}".format(sum(d_nonccs.itervalues())))
    print("Avg. CCS readlength: {0}".format(sum(r_ccs.itervalues())/a))
    print("Avg. non-CCS subread readlength: {0}".format(sum(r_nonccs.itervalues())/len(r_nonccs)))
    
def summarize_chimera(ccs_prefix='ccs_reads.53Aseen_trimmed_changeid.fa', nonccs_prefix='nonccs_subreads.53Aseen_trimmed_changeid.fa'):
    """
    1. % of artificial chimeras in CCS
    2. % of artificial chimeras in subreads
    3. % of ....combined
    """
    d_ccs_is = get_zmw(ccs_prefix + '.is_chimera.fa')
    d_ccs_not = get_zmw(ccs_prefix + '.non_chimera.fa')
    d_nonccs_is = get_zmw(nonccs_prefix + '.is_chimera.fa')
    d_nonccs_not = get_zmw(nonccs_prefix + '.non_chimera.fa')    
    
    a = len(d_ccs_is)
    b = a + len(d_ccs_not)
    print("% of artificial chimeras in CCS: {0}/{1} ({2:.1f}%)").format(a, b, a*100./b)
    a = len(d_nonccs_is)
    b = a + len(d_nonccs_not)
    print("% of artificial chimeras in non-CCS subreads: {0}/{1} ({2:.1f}%)").format(a, b, a*100./b)
    a = len(d_ccs_is) + len(d_nonccs_is)
    b = a + len(d_ccs_not) + len(d_nonccs_not)
    print("% of artificial chimeras: {0}/{1} ({2:.1f}%)").format(a, b, a*100./b)
    
def summarize_gmap(ccs_prefix='ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa', nonccs_prefix='nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa'):    
    ccs_unmapped, ccs_obs2exp, (ccs_chimera_missed,ccs_chimera_real,ccs_chimera_rate), ccs_avg_coverage, ccs_zmw_count = eval_accuracy_by_gmap.run_gmap(ccs_prefix, None, ccs_prefix+'.gff')
    nonccs_unmapped, nonccs_obs2exp, (nonccs_chimera_missed,nonccs_chimera_real,nonccs_chimera_rate), nonccs_avg_coverage, nonccs_zmw_count = eval_accuracy_by_gmap.run_gmap(nonccs_prefix, None, nonccs_prefix+'.gff')
        
    ccs_avg_obs = ccs_obs2exp['ObsAccuracy'].mean()
    nonccs_avg_obs = nonccs_obs2exp['ObsAccuracy'].mean()
    
    print "----- CCS only ------"
    print "Total number of ZMWs:", ccs_zmw_count
    print "Total number of unmapped: {0} ({1:.1f}%)".format(ccs_unmapped, ccs_unmapped*100./ccs_zmw_count)
    print "Avg. coverage: {0:.1f}%".format(ccs_avg_coverage)
    print "Avg. observed accuracy: {0:.1f}%".format(ccs_avg_obs)
    print "Chimera rate: {0:.1f}%".format(ccs_chimera_rate*100.)
    print "Chimera from missed adapter: {0:.1f}%".format(ccs_chimera_missed*100./(ccs_chimera_missed+ccs_chimera_real))
    print "----- non-CCS subread only ------"
    print "Total number of ZMWs:", nonccs_zmw_count
    print "Total number of unmapped: {0} ({1:.1f}%)".format(nonccs_unmapped, nonccs_unmapped*100./nonccs_zmw_count)
    print "Avg. coverage: {0:.1f}%".format(nonccs_avg_coverage)
    print "Avg. observed accuracy: {0:.1f}%".format(nonccs_avg_obs)
    print "Chimera rate: {0:.1f}%".format(nonccs_chimera_rate*100.)
    print "Chimera from missed adapter: {0:.1f}%".format(nonccs_chimera_missed*100./(nonccs_chimera_missed+nonccs_chimera_real))
    print "----- Combined -------"
    a = ccs_unmapped + nonccs_unmapped
    b = ccs_zmw_count + nonccs_zmw_count
    print "Total number of ZMWs:", b
    print "Total number of unmapped: {0} ({1:.1f}%)".format(a, a*100./b)
    avg_coverage = (ccs_avg_coverage*ccs_zmw_count + nonccs_avg_coverage*nonccs_zmw_count)*1./(ccs_zmw_count + nonccs_zmw_count)
    print "Avg. coverage: {0:.1f}%".format(avg_coverage)
    avg_obs = (ccs_avg_obs*ccs_zmw_count + nonccs_avg_obs*nonccs_zmw_count)*1./(ccs_zmw_count + nonccs_zmw_count)
    print "Avg. observed accuracy: {0:.1f}%".format(avg_obs)
    chimera_rate = (ccs_chimera_rate*ccs_zmw_count + nonccs_chimera_rate*nonccs_zmw_count)*1./(ccs_zmw_count + nonccs_zmw_count)
    print "Chimera rate: {0:.1f}%".format(chimera_rate*100.)

def parse_blasr(sam_filename, ref_fasta_filename):
    """
    Return dict of ZMW --> best r by maximizing sCov
    """
    hit = {}
    for r in miscBio.SAMReader(sam_filename, True, ref_fasta_filename):
        zmw = r.qID[:r.qID.find('/', r.qID.find('/')+1)]
        if zmw not in hit or r.sCoverage >= hit[zmw].sCoverage:
            hit[zmw] = r
    return hit

def tally_ref_hits(hit_by_zmw, ref_tally, min_sCov, min_qCov):
    for r in hit_by_zmw.itervalues():
        if r.sCoverage >= min_sCov and r.qCoverage >= min_qCov:
            ref_tally[r.sID] += 1

def summarize_blasr(ccs_prefix='ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa',nonccs_prefix='nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa', ref_fasta_filename='/home/UNIXHOME/etseng/share/gencode/gencode.v15.pc_transcripts.non_redundant_good.fa.nonredundant.fasta', ref_size=None):
    """
    1. Number of aligned ZMWs / refs
    2. Number of well-aligned (90% qCov, sCov) ZMWs / refs
    3. Abundance distribution of well-aligned refs
    """
    hit_ccs = eval_refmap_by_blasr.parse_blasr(ccs_prefix+'.blasr.sam', ref_fasta_filename)
    hit_nonccs = eval_refmap_by_blasr.parse_blasr(nonccs_prefix+'.blasr.sam', ref_fasta_filename)

    tally0 = defaultdict(lambda: 0)
    tally90 = defaultdict(lambda: 0)

    eval_refmap_by_blasr.tally_ref_hits(hit_ccs, tally0, 0., 0.)
    eval_refmap_by_blasr.tally_ref_hits(hit_nonccs, tally0, 0., 0.)
    eval_refmap_by_blasr.tally_ref_hits(hit_ccs, tally90, 0.9, 0.9)
    eval_refmap_by_blasr.tally_ref_hits(hit_nonccs, tally90, 0.9, 0.9)

    ref_count = int(os.popen("grep -c \">\" " + ref_fasta_filename).read().strip())
    zmw_count = len(get_zmw(ccs_prefix)) + len(get_zmw(nonccs_prefix))

    a = len(tally0)
    b = ref_count
    c = len(hit_ccs) + len(hit_nonccs)
    d = zmw_count
    print("Number of aligned refs: {0}/{1} ({2:.1f}%)".format(a, b, 100.*a/b))
    print("Number of aligned ZMWs: {0}/{1} ({2:.1f}%)".format(c, d, 100.*c/d))
    a = len(tally90)
    c = sum(r.sCoverage>=.9 and r.qCoverage>=.9 for r in hit_ccs.itervalues()) +\
        sum(r.sCoverage>=.9 and r.qCoverage>=.9 for r in hit_nonccs.itervalues())
    print("Number of well-aligned refs: {0}/{1} ({2:.1f}%)".format(a, b, 100.*a/b))
    print("Number of well-aligned ZMWs: {0}/{1} ({2:.1f}%)".format(c, d, 100.*c/d))

    # draw plots
    hit = hit_ccs
    hit.update(hit_nonccs)

    eval_refmap_by_blasr.draw_2dhist(hit, 'ccs_n_nonccs.sLen_vs_qLen', feat_func=lambda x: (x.sLen, x.qLen), filter_func=lambda x: True, xlab='Reference Length', ylab='Query length')
    eval_refmap_by_blasr.draw_2dhist(hit, 'ccs_n_nonccs.sLen_vs_qLen_well', feat_func=lambda x: (x.sLen, x.qLen), filter_func=lambda x: x.qCoverage>=.9 and x.sCoverage>=.9, xlab='Reference Length', ylab='Query length')
    eval_refmap_by_blasr.draw_2dhist(hit, 'ccs_n_nonccs.sCov_vs_qCov', feat_func=lambda x: (x.sCoverage, x.qCoverage), filter_func=lambda x: True, xlab='Reference Coverage', ylab='Query Coverage')
    eval_refmap_by_blasr.draw_2dhist(hit, 'ccs_n_nonccs.sLen_vs_sCov', feat_func=lambda x: (x.sLen, x.sCoverage), filter_func=lambda x: True, xlab='Reference Length', ylab='Reference Coverage')

    if ref_size is None:
        eval_refmap_by_blasr.draw_coverage(hit, 'ccs_n_nonccs.ref_coverage', filter_func=lambda x: x.qCoverage>=.9, title="5'-3' reference coverage (qCov>=90%)")
    else:
        a, b = ref_size
        eval_refmap_by_blasr.draw_coverage(hit, 'ccs_n_nonccs.ref_coverage', filter_func=lambda x: x.qCoverage>=.9 and a<=x.sLen<=b, title="5'-3' reference coverage (qCov>=90%), ref length {0}-{1} bp".format(a,b))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--ref_fasta_filename", default='/home/UNIXHOME/etseng/share/gencode/gencode.v15.pc_transcripts.non_redundant_good.fa.nonredundant.fasta')
    parser.add_argument("--ref_size", default=None)
    args = parser.parse_args()


    summarize_CCS()
    summarize_chimera()
    summarize_gmap()
    summarize_blasr(ref_fasta_filename=args.ref_fasta_filename, ref_size=args.ref_size)
    