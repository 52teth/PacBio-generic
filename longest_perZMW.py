import os, sys
import numpy as np
from cPickle import *

"""
Suite of functions that is based on longest subread per ZMW
"""

def iter_longest_perZMW(data):
    """
    alns is the pickle output from alignQC
    could also be inserts, in which case it'd be ALL ZMWs
    if it's alns, then it's just subreads that aligned

    for each ZMW (key: MovieID+HoleNumber), select the representative
    to be the longest subread. 
    """
    if type(data) is dict: # is inserts
        # need to change rStart,rEnd to iStart,iEnd
        for movie, inss in data.iteritems():
            names = inss.dtype.names
            assert names[1] == 'rStart' and names[2] == 'rEnd'
            inss.dtype.names = ('HoleNumber','iStart','iEnd','IsFullPass','IsHQTrimmed','IsLongest','IsAT')
            holes = np.unique(inss['HoleNumber'])
            for hn in holes:
                yield movie, inss[(inss['HoleNumber']==hn)&(inss['IsLongest'])]            
    else: # is alns
        alns_whole = data
        movies = np.unique(alns_whole['MovieID'])
        for movie in movies:
            alns = alns_whole[alns_whole['MovieID']==movie]
            holes = np.unique(alns['HoleNumber'])
            for hn in holes:
                yield movie, alns[(alns['HoleNumber']==hn)&(alns['IsLongest'])]            

def report_gene_match(alns_whole, output_filename):
    """
    Report gene match per longest subread per ZMW
    """
    f = open(output_filename, 'w')
    f.write("MovieID\tHoleNumber\tisAt\tiCoverage\ttName\ttStart\ttEnd\n")
    for a in iter_longest_perZMW(alns_whole):
        if len(a) == 0: continue
        cov = (a['rEnd']-a['rStart'])*1./(a['iEnd']-a['iStart'])
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(\
                a['MovieID'][0], a['HoleNumber'][0], a['IsAT'][0],\
                cov[0], a['RefID'][0], a['tStart'][0], a['tEnd'][0]))
    f.close()

def extract_longest_subread_with_bad_coverage(alns_whole, movieID_dict, subreads_filename, coverage_cutoff=0.9):
    to_use = []
    for a in iter_longest_perZMW(alns_whole):
        if len(a) == 0: continue
        cov = (a['rEnd']-a['rStart'])*1./(a['iEnd']-a['iStart'])
        if cov[0] < coverage_cutoff:
            seqid = "{movie}/{hole}/{s}_{e}".format(movie=movieID_dict[a['MovieID'][0]],\
                    hole=a['HoleNumber'][0], s=a['iStart'][0], e=a['iEnd'][0])
            to_use.append(seqid)
            #print seqid
            #raw_input()

    from Bio import SeqIO
    for r in SeqIO.parse(open(subreads_filename), 'fasta'):
        if r.id in to_use:
            to_use.remove(r.id)
            print(">{0}\n{1}".format(r.id, r.seq))

    print >> sys.stderr, to_use
    
def extract_longest_subread(alns_whole, subreads_filename, output_filename):
    """
    alns_whole --- aln or ins output from alignQC
    movieID_dict --- movieID number (ex: 1) to movieName
    subreads_filename --- filtered_subreads.fasta
    """
    to_use = []
    for movieID,a in iter_longest_perZMW(alns_whole):
        if len(a) == 0: continue
        seqid = "{movie}/{hole}/{s}_{e}".format(movie=movieID,\
                    hole=a['HoleNumber'][0], s=a['iStart'][0], e=a['iEnd'][0])
        to_use.append(seqid)
        
    from Bio import SeqIO
    f = open(output_filename, 'w')
    for r in SeqIO.parse(open(subreads_filename), 'fasta'):
        if r.id in to_use:
            to_use.remove(r.id)
            f.write(">{0}\n{1}\n".format(r.id, r.seq))
    f.close()

    print >> sys.stderr, to_use    

if __name__ == "__main__":
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("--pickle", required=True, help="Pickle filename (from alignQC)")
    parser.add_argument("--output", required=True, help="Output filename")
    parser.add_argument("--subreads", help="Filtered subreads (optional; needed by extract longest subread)")
    parser.add_argument("--action", required=True, choices=['extract_longest'])

    
    args = parser.parse_args()
    
    if args.action == 'extract_longest':
        assert args.subreads is not None

    with open(args.pickle) as f:
        ins, aln = load(f)

    movieID_dict = {
            1:       'm121103_033957_42175_c100419402550000001523040112191250_s1_p0',
            2:       'm121103_055923_42175_c100419402550000001523040112191251_s1_p0',
            3:       'm121103_130510_42175_c100419402550000001523040112191254_s1_p0',
            4:       'm121103_152831_42175_c100419402550000001523040112191255_s1_p0'}
    subreads_filename = '/home/UNIXHOME/etseng/projects/swiss_neurexin/smrtpipe/library1/data/filtered_subreads.fasta'

    #extract_longest_subread_with_bad_coverage(aln, movieID_dict, subreads_filename, 0.8)
    if args.action == 'extract_longest':
        extract_longest_subread(ins, args.subreads, args.output)