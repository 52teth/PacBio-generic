import os, sys
import random
import bisect
from cPickle import *
from collections import defaultdict
import numpy as np

def rarefaction_sample(acc_counts, size, iters=1):
    N = acc_counts[-1]
    assert size <= N
    result = []
    for iter in xrange(iters):
        seen = set()
        for s in random.sample(range(N), size):
            seen.add(bisect.bisect_left(acc_counts, s+1))
        result.append(len(seen))
    return result

def main(pickle_filename):
    """
    input should be a pickle file from running alignQC
    must have key 'alns' and 'refDict'
    """
    with open(pickle_filename) as f:
        p = load(f)
        aln = p['alns']

    # use only IsAT & refCoverage >= 80% & per-ZMW
    good = aln[((aln['rEnd']-aln['rStart'])>=0.8*(aln['iEnd']-aln['iStart']))&((aln['tEnd']-aln['tStart'])>=0.8*aln['RefLength'])&aln['IsAT']&aln['IsFullPass']]
    seen_moviehole = set()
    counts = defaultdict(lambda: 0) # RefID --> counts
    for x in good:
        moviehole = x['MovieID'],x['HoleNumber']
        if moviehole in seen_moviehole: continue
        seen_moviehole.add(moviehole)
        counts[x['RefID']] += 1

    # change counts into accumulative form
    counts = counts.values()
    acc = [counts[0]]
    for c in counts[1:]: acc.append(acc[-1]+c)

    # do rarefaction
    N = acc[-1]
    for size in xrange(0, N, 100):
        sampled = rarefaction_sample(acc,size,iters=50)
        print size, np.mean(sampled), np.std(sampled)

if __name__ == "__main__":
    main(sys.argv[1])
