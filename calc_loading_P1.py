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
        f.write("Name,allZmw,seqZmw,P1\n")
        for file in os.popen("ls {0}/*/Analysis_Results/*.bas.h5".format(run_dir)).read().strip().split('\n'):
            bas = BasH5Reader(file)
            name = os.path.basename(file)
            a = len(bas.sequencingZmws)
            b = len(bas.allSequencingZmws)
            f.write("{0},{1},{2},{3:.2f}\n".format(bas.movieName, b, a, 100.*a/b))
            result[bas.movieName] = (a, b)
        f.close()
    return result

if __name__ == "__main__":
    get_P1(os.path.abspath(sys.argv[1]))

