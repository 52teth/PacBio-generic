#!/usr/bin/env python
import os, sys, fnmatch
from csv import DictReader
from pbcore.io import BasH5Reader

def get_P1(run_dir, force_recalc=False):
    result = {} # movieName --> (# of sequencing ZMWs, # of ZMWs)
    
    filename = os.path.join(run_dir, 'loading_P1.txt')
    if os.path.exists(filename) and not force_recalc:
        with open(filename) as f:
            for r in DictReader(f, delimiter=','):
                result[r['Name']] = (int(r['seqZmw']), int(r['allZmw']))                
    else:
        f = open(filename, 'w')
        f.write("Name,allZmw,seqZmw,P1,avgLen\n")
        for file in os.popen("ls {0}/*/Analysis_Results/*.bax.h5".format(run_dir)).read().strip().split('\n'):
            bas = BasH5Reader(file)
            name = os.path.basename(file)
            a = len(bas.sequencingZmws)
            b = len(bas.allSequencingZmws)
            lens = [len(bas[x].read()) for x in bas.sequencingZmws]
            f.write("{0},{1},{2},{3:.2f},{4:.0f}\n".format(bas.movieName, b, a, 100.*a/b, sum(lens)*1./a))
            if bas.movieName not in result: result[bas.movieName] = (a, b)
            else: result[bas.movieName] = (a+result[bas.movieName][0], b+result[bas.movieName][1])
        f.close()
    return result

if __name__ == "__main__":
    get_P1(os.path.abspath(sys.argv[1]))

