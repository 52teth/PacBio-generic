#!/usr/bin/env python

import os, sys
#sys.path.insert(0, '/home/UNIXHOME/etseng/.VENV/lib/python2.7/site-packages')
#print >> sys.stderr, "IMPORTED VENV for matplotlib 1.2.0!!!"
import numpy as n
from re import sub
import argparse
import cPickle
import pdb
import itertools

from pbcore.io.BasH5IO import BasH5
from pbcore.io.cmph5 import factory, alignmentPairMap

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter, LinearLocator

from scipy.stats import mode, mstats

from csv import DictReader
from collections import defaultdict
from bx.intervals.intersection import IntervalTree

def getMaxCoveredRefMatch(alns, refLengthDict):
    """
    alns is a slice of /AlnInfo
    (could be all alignments within a subread or ZMW)
    
    Return just the one that covers the max proportion of the ref (not the read)  
    """
    return max(alns, key=lambda x: (x['tEnd']-x['tStart'])*1./refLengthDict[x['RefGroupID']])

def group_by_inserts(alns, inss):
    """
    alns is a slice of /AlnInfo
    inss is the corressponding slice from getInsertsFromBasH5
    
    return a list of ins_index --> list of aln_index within this insert
    """
    grouped = [[] for i in xrange(len(inss))]
    for aln_index, aln in enumerate(alns):
        found = False
        for ins_index, ins in enumerate(inss):
            if ins['rStart'] <= aln['rStart'] <= aln['rEnd'] <= ins['rEnd']:
                grouped[ins_index].append(aln_index)
                found = True
                break
        if not found:
            raise Exception, "alignment not found within inserts!!! Not OK!!!"
    
    # if doing something weird like partial alignment, comment the raise Exception above and use this!
#    i = 0
#    while i < len(grouped):
#        if len(grouped[i]) == 0:
#            grouped.pop(i)
#        else: i += 1

    return grouped

def getPrimerInfo(filename):
    """
    Read .primer_info.txt from running parse_seq_clean.py
    <ID> <5seen> <polyAseen> <3seen>
    
    Returns: dict of (movie, hn) --> IntervalTree with dict record of {'ID', '5seen', 'polyAseen', '3seen'}
    """
    d = defaultdict(lambda: IntervalTree())
    with open(filename) as f:
        for r in DictReader(f, delimiter='\t'):
            movie, hn, s_e = r['ID'].split('/')
            s, e = map(int, s_e.split('_'))
            d[(movie, hn)].insert(s, e, r)
    return d

# OBSOLETE! Use PrimerInfo
#def getPolyAT(filename):
#    """
#    Read polyAT report file (from seqclean)
#    Movie   HoleNumber      rStart  rEnd
#    
#    Returns a dict of (movieName, holeNumber) --> list of (start, end) of polyAT
#    """
#    from csv import DictReader
#    d = defaultdict(lambda: [])
#    for r in DictReader(open(filename), delimiter='\t'):
#        d[(r['Movie'],r['HoleNumber'])].append((int(r['rStart']), int(r['rEnd'])))
#    return d
        
def hasPolyAT(d, movie, hn, start, end):
    """
    Best used with functools.partial to create a boolean function for checking
    whether an insert is in the dictionary
    
    NOTE: now modified to work with primer_info.txt
    NOTE: returns True iff 5seen=='1' and 3seen=='1' !!! polyA actually ignored.
    """
    k = movie, hn
    if k not in d:
        return False
    
    item = d[k].find(start, end)
    if len(item) == 0:
        return False
    elif len(item) == 1:
        rec = item[0]
        assert rec['5seen'] in ('0', '1') and rec['3seen'] in ('0', '1')
        return rec['5seen'] == '1' and rec['3seen'] == '1'
    else:
        print >> sys.stderr, "Impossible to have {0}/{1}/{2}_{3}!! I probably hacked something".format(movie, hn, start, end)
  

def getInsertsFromBasH5(bash5FN, func_is_polyAT):
    """
    Reads through a single .bas.h5 and return an array of
    
    HoleNumber, rStart, rEnd, IsFullPass, IsHQTrimmed, IsLongest, IsAT
    
    NOTE: IsLongest not always equate to IsFullPass esp. with StageStart
    NOTE: IsAT is actually set according to 5seen & 3seen. Check function definition!!
    """
    bash5 = BasH5(bash5FN)
    rgnTable = bash5.rgnTable
    data = []    

    if 'Adapter' in bash5.rgnTable.rgnDS.attrs.get('RegionTypes',[]):
        for hn in bash5.getSequencingZMWs():
            # get indices first
            adapters = rgnTable.getAdapterRegionForZMW(hn)
            hqregion = rgnTable.getHQRegionForZMW(hn)
            inserts = rgnTable.getInsertRegionForZMW(hn)

            adapterStarts = [n.min(x) for x in adapters]
            adapterEnds = [n.max(x) for x in adapters]
            hqStart = hqregion[0][0]
            hqEnd = hqregion[0][1]
            
            hqinserts = filter(lambda ins: hqStart<=ins[0]<hqEnd or hqStart<=ins[1]<hqEnd, inserts)
            if len(hqinserts) == 0:
                continue

            for i in hqinserts:
                iStart, iEnd = list(i)

                 # check if insert is full pass
                isFullPass = iStart in adapterEnds and iEnd in adapterStarts                
                isLongest = False # for now set everything to false, re-set at the end                
                isHQTrimmed = False
                # trim by HQRegion
                if iStart < hqStart:
                    iStart = hqStart
                    isHQTrimmed = True
                if iEnd > hqEnd:
                    iEnd = hqEnd
                    isHQTrimmed = True

                # this has to be after HQ trim to be correct                    
                isAT = func_is_polyAT(str(hn), iStart, iEnd)

                insert = (hn, iStart, iEnd, isFullPass, isHQTrimmed, isLongest, isAT)
                #print "Adding", insert
                data.append(insert)

    data = n.array(data, dtype=[('HoleNumber', '<i4'), ('rStart', '<i4'), ('rEnd', '<i4'), ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool), ('IsAT', n.bool)])
    # for each ZMW, set the IsLongest label                
    for hn in n.unique(data['HoleNumber']):
        x = data[data['HoleNumber']==hn]
        max_len = max(x['rEnd'] - x['rStart'])
        p = data[(data['HoleNumber']==hn)&(data['rEnd']-data['rStart']==max_len)]
        p['IsLongest'] = True
        data[(data['HoleNumber']==hn)&(data['rEnd']-data['rStart']==max_len)] = p
    return data
        
def getInsertsFromFofn(inputFOFN, primer_match_filename):
    """
    primer_match_filename should be .primer_info.txt
    
    Returns: dict of movieName --> array from reading .bas.h5
    """
    import functools    
    primer_match_dict = getPrimerInfo(primer_match_filename)
    inserts = {}
    
    with open(inputFOFN) as f:
        for line in f:
            filename = line.strip()
            print "Reading", filename
            movieName = sub('.pls.h5|.bas.h5', '', os.path.basename(filename))
            inserts[movieName] = getInsertsFromBasH5(filename, functools.partial(hasPolyAT, primer_match_dict, movieName))

    return inserts

def getReferenceLengths(cmph5):
    return cmph5['/RefInfo'].asRecArray()['Length']


def getAlignedLengthRatios(cmph5, inserts):
    """
    Expects that BLASR could've been run with -bestN > 1.
    Given multiple ref hits, reports just the one with the *max ref coverage*
    """
    aIdx = cmph5['/AlnInfo'].asRecArray()
    refLengthDict = cmph5['/RefInfo'].asDict("ID", "Length", cache=True)
    refIdDict = cmph5['/RefGroup'].asDict("ID", "RefInfoID", cache=True)

    results = []

    movieDict = dict(zip(cmph5['/MovieInfo/Name'], cmph5['/MovieInfo/ID']))
    for movieName in inserts:
        if movieName in movieDict:
            movieID = movieDict[movieName]
            ins = inserts[movieName]
            sl = aIdx[aIdx['MovieID'] == movieID]

            # get aligned reads' hole numbers
            for hn in n.unique(sl['HoleNumber']):

                # get inserts and alignments
                inss = ins[ins['HoleNumber'] == hn]
                alns = sl[sl['HoleNumber'] == hn]
                
                if len(inss) == 0:
                    continue
                
                for ins_index, aln_indices in enumerate(group_by_inserts(alns, inss)):
                    if len(aln_indices) == 0:
                        continue
                    i = inss[ins_index]
                    a = getMaxCoveredRefMatch(alns[aln_indices], refLengthDict)
                   # pdb.set_trace()
                    rStart = a['rStart']
                    rEnd = a['rEnd']
                    tStart = a['tStart']
                    tEnd = a['tEnd']
                    refGroupID = a['RefGroupID']
                    refLength = refLengthDict[refIdDict[refGroupID]]
                    iStart = i['rStart']
                    iEnd = i['rEnd']
                    iIsFullPass = i['IsFullPass']
                    iIsHQTrimmed = i['IsHQTrimmed']
                    iIsLongest = i['IsLongest']
                    iIsAT = i['IsAT']
                    refStrand = a['RCRefStrand']

                    results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, tStart, tEnd, iStart, iEnd, refLength, iIsFullPass, iIsHQTrimmed, iIsLongest, iIsAT, refStrand))

    return n.array(results, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),
                                   ('rStart', '<i4'), ('rEnd', '<i4'),
                                   ('tStart', '<i4'), ('tEnd', '<i4'),
                                   ('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'),
                                   ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool),
                                   ('IsAT', n.bool), ('refStrand', '<i4')])
                                        

def makeFractionSubreadHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Subread in Alignment Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    Longest = alnRatios[alnRatios['IsLongest']]
    AT = alnRatios[alnRatios['IsAT']]

    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):
        alnLength = l['rEnd'] - l['rStart']
        srLength = l['iEnd'] - l['iStart']
        if len(srLength) == 0:
            continue
        alnSubRatio = alnLength.astype(float) / srLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper left')
    ax.set_xlabel("Alignment length/Subread Length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeSubreadRLHistogram(alnRatios, outfile, format, quantile=None, ylim=None):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Subread RL Histogram")
    
    if ylim is not None:
        plt.ylim(ylim)

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    Longest = alnRatios[alnRatios['IsLongest']]
    AT = alnRatios[alnRatios['IsAT']]
    
    num_bins_fordraw = max(Longest['iEnd']-Longest['iStart']) / 50. # bins by 50bp
    
    max_y = 0
    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):
        srLength = l['iEnd'] - l['iStart']
        if len(srLength) == 0:
            continue
        if quantile is not None:
            srLength = srLength[srLength < mstats.mquantiles(srLength, [quantile])[0]]
        num, ignore1, ignore2 = ax.hist(srLength, bins=100, histtype='step', label=label)
        max_y = max(max(num), max_y)
        
    ax.set_ylim(0, max_y*1.1)
    ax.legend(loc='upper right', ncol=1)
    ax.set_xlabel("Subread Length")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeFractionReferenceHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Reference in Alignment Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    Longest = alnRatios[alnRatios['IsLongest']]
    AT = alnRatios[alnRatios['IsAT']]
    
    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):    
        alnLength = l['tEnd'] - l['tStart']
        if len(alnLength) == 0:
            continue
        refLength = l['RefLength']
        alnSubRatio = alnLength.astype(float) / refLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper right')
    ax.set_xlabel("Alignment Length/Reference Length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeReferenceRLHistogram(alnRatios, refLengths, outfile, format, quantile=None):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned References RL Density Plot")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]

    max_y = 0
    for l, label in zip((alnRatios, HQRegion, fullPass), ("Aln from All Subreads", "Aln from HQRegion Subreads", "Aln from HQRegion Full-Pass Subreads")):
        alnRefLength = l['RefLength']
        if len(alnRefLength) == 0:
            continue

        if not quantile == None:
            alnRefLength = alnRefLength[alnRefLength < mstats.mquantiles(alnRefLength, [quantile])[0]]

        num, bins, patches = ax.hist(alnRefLength, bins=100, histtype='step', label=label, normed=True)

        if n.max(num) > max_y:
            max_y = n.max(num)

    if not quantile == None:
        refLengths = refLengths[refLengths < mstats.mquantiles(refLengths, [quantile])[0]]

    num, bins, patches = ax.hist(refLengths, bins=100, histtype='step', label="All References", normed=True)

    if n.max(num) > max_y:
        max_y = n.max(num)

    ax.set_ylim(0, max_y * 1.1)
    ax.legend(loc='upper center', prop={'size': 'small'})
    ax.set_xlabel("Reference Length")
    ax.set_ylabel("Density")
    fig.savefig(outfile, format=format)

def makeAlignmentPercentileDistribution(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refLengthRange=None, per_gene=False, refStrandDict=None):
    """
    This is the 5'-3' reference coverage plot
    
    refLengthRange --- plot only for genes with refLength in range
    per_gene --- plot coverage using only the max coverage per gene 
    """
    fig = plt.figure(dpi=300, figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    #ax1.set_title("Normalized Reference Start/End Positions")
    title = "Reference coverage ({0} only, qStart<100)".format(label)
    if per_gene:
        title += ",per_gene"
    ax1.set_title(title)   

    ax2 = fig.add_subplot(212)
    #ax2.set_title("Normalized Reference Coverage")
    if refLengthRange is None:
        ax2.set_title("Reference Coverage")
    else:
        ax2.set_title("Reference Coverage (size:{0}-{1})".format(refLengthRange[0], refLengthRange[1]))
    
    if refLengthRange is None:
        fullPass = alnRatios[alnRatios[alnKey]&(alnRatios['rStart']-alnRatios['iStart']<100)]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(alnRatios['rStart']-alnRatios['iStart']<100)&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]
        
    if len(fullPass) == 0:
        print >> sys.stderr, "Not drawing reference coverage for {0}!!".format(label)
        return
        
    if per_gene:
        refLength = []
        startPercentile = []
        endPercentile = []
        for refID in n.unique(fullPass['RefID']):
            max_cov = 0
            reflen = 0
            for x in fullPass[fullPass['RefID']==refID]:
                reflen = x['RefLength']*1.
                cov = x['tEnd'] - x['tStart']
                if cov > max_cov:                    
                    max_cov = cov
                    s = x['tStart'].astype(float) / reflen
                    e = x['tEnd'].astype(float) / reflen
            assert max_cov > 0
            refLength.append(reflen)
            if refStrandDict is None or refStrandDict[refID]=='+':
                startPercentile.append(s)
                endPercentile.append(e)
            else:
                startPercentile.append(1-e)
                endPercentile.append(1-s)
                            
    else:    
        if refStrandDict is None:            
            refLength = fullPass['RefLength']
            startPercentile = fullPass['tStart'].astype(float) / refLength.astype(float)
            endPercentile = fullPass['tEnd'].astype(float) / refLength.astype(float)
        else:
            startPercentile = []
            endPercentile = []
            for x in fullPass:
                reflen = x['RefLength']*1.
                s = x['tStart'].astype(float) / reflen
                e = x['tEnd'].astype(float) / reflen
                if refStrandDict[x['RefID']] == '+':
                    startPercentile.append(s)
                    endPercentile.append(e)
                else:
                    startPercentile.append(1-e)
                    endPercentile.append(1-s)
                    

    max_y = 0
    for l, label, color in itertools.izip([startPercentile, endPercentile], ['Alignment Start', 'Alignment End'], ('#FF33CC', '#0000CC')):
        num, bins, patches = ax1.hist(l, bins=101, histtype='stepfilled', label=label, normed=True, color=color)
        max_y = max(max_y, n.max(num))


    acc = n.zeros(101)
    for a,b in itertools.izip(startPercentile,endPercentile): acc[int(a*100):int(b*100)+1] += 1
    ax2.fill_between(n.arange(0,1.01,0.01), acc, facecolor='#3399FF', alpha=.7)
    
    ax2.set_xlabel("Normalized Reference Sequence Length")
    ax2.set_ylabel("Count")
                                 
    ax1.set_ylim(0, max_y * 1.1)
    #ax1.set_xlim(-0.01, 1.01)
    ax1.legend(loc='upper center', prop={'size': 'small'})
    
    fig.savefig(outfile, format=format)
    
    
def makeCoverage_by_RefLength(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass'):    
    title = "Reference Length vs Normalized Reference Coverage ({0} only)".format(label)
    fullPass = alnRatios[alnRatios[alnKey]]
    refLength = fullPass['RefLength']
    startPercentile = fullPass['tStart'].astype(float) / refLength.astype(float)
    endPercentile = fullPass['tEnd'].astype(float) / refLength.astype(float)
    
    X = n.round(refLength, decimals=-2)
    Y = n.round(endPercentile - startPercentile, decimals=2)
    
    _makeHexbinHist(X, Y, "Reference Length", "Normalized Reference Coverage", title, outfile, format, quantile=0.99)

def makeStartPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refStrandDict=None):
    title = label + " Subreads vs Reference alignment start"
    fullPass = alnRatios[alnRatios[alnKey]]
    
    if refStrandDict is None:
        X = fullPass['tStart'].astype(float)
        Y = (fullPass['rStart'] - fullPass['iStart']).astype(float)
    else:
        X = []
        Y = (fullPass['rStart'] - fullPass['iStart']).astype(float)
        for x in fullPass:
            if refStrandDict[x['RefID']] == '+': X.append(x['tStart'].astype(float))
            else: X.append((x['RefLength']-x['tEnd']).astype(float))
        X = n.array(X)
                    
    _makeHexbinHist(X, Y, "Reference alignment start", "Subread alignment start", title, outfile, format, quantile=0.99)

def makeEndPositionVS(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refStrandDict=None):
    title = label + " Subreads vs Reference alignment distance from end"
    fullPass = alnRatios[alnRatios[alnKey]]
    if refStrandDict is None:
        X = (fullPass['RefLength'] - fullPass['tEnd']).astype(float)
        Y = (fullPass['iEnd'] - fullPass['rEnd']).astype(float)
    else:
        X = []
        Y = (fullPass['iEnd'] - fullPass['rEnd']).astype(float)
        for x in fullPass:
            if refStrandDict[x['RefID']] == '+': X.append((x['RefLength']-x['tEnd']).astype(float))
            else: X.append(x['tStart'].astype(float))
        X = n.array(X)
    
    _makeHexbinHist(X, Y, "Distance from ref end", "Distance from subread end", title, outfile, format, quantile=0.99)
    
#def makeRefAbundance(alnRatios, outfile, format):
#    fig = plt.figure(dpi=300, figsize=(10, 6))
#    ax = fig.add_subplot(111)
#    ax.set_title("Abundance of aligned references (>= 80% ref aligned)")   
#    
#    fullPass = alnRatios[alnRatios['IsFullPass']]
#    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
#    Longest = alnRatios[alnRatios['IsLongest']]
#    AT = alnRatios[alnRatios['IsAT']]
#    
#    max_y = 0
#    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):
#        count = {} # refID --> hits
#        for x in l:
#            cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
#            if cov < .8:
#                continue
#            if x['RefID'] not in count: count[x['RefID']] = 1
#            else: count[x['RefID']] += 1
#        x_max = max(count.itervalues())
#        num, ignore1, ignore2 = ax.hist(count.values(), bins=x_max/2, histtype='step', label=label)
#        max_y = max(max(num), max_y)
#        
#    ax.set_ylim(0, max_y*1.1)
#    ax.legend(loc='upper right', ncol=1)
#    ax.set_xlabel("Number of subread hits")
#    ax.set_ylabel("Count")
#    fig.savefig(outfile, format=format) 

    
def makeRefLengthVSAbundance(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', covThreshold=.8):
    title = "Reference Length vs Subread hits ({0} only, >= {1}% ref aligned)".format(label, n.round(covThreshold*100))
    fullPass = alnRatios[alnRatios[alnKey]]

    count = {} # refID --> hits
    reflen_dict = {} # refID --> refLength
    for x in fullPass:
        cov = (x['tEnd'] - x['tStart'])*1. / x['RefLength']
        if cov < covThreshold: 
            continue
        if x['RefID'] not in count:
            count[x['RefID']] = 1
            reflen_dict[x['RefID']] = n.round(x['RefLength'], decimals=-2) # round refLength by 100bp
        else: count[x['RefID']] += 1
        
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title(title)
    
    from collections import defaultdict
    d = defaultdict(lambda: 0) # (ref length, abundance) --> frequency
    for refID, hits in count.iteritems():
        d[(reflen_dict[refID], hits)] += 1
        
    X, Y, Z = [], [], []
    for (x,y),z in d.iteritems():
        X.append(x)
        Y.append(y)
        Z.append(z*20)
        #print >> sys.stderr, x, y, z
    
    ax.scatter(X, Y, Z, alpha=0.7, marker='o', facecolor='#3399FF')
    ax.set_xlabel("Reference Length")
    ax.set_ylabel("Number of subread hits")
    ax.set_xlim(0, max(X)+100)
    ax.set_ylim(0, max(Y)*1.1)
    fig.savefig(outfile, format=format)
        
def makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', title='Full-Pass'):
    title = title + " Subread Length vs Reference Length"
    x_label = "Reference Length"
    y_label = "Subread Length"

    fullPass = alnRatios[alnRatios[alnKey]]
    subreadLength = fullPass['iEnd'] - fullPass['iStart']
    refLength = fullPass['RefLength']
    
    _makeHexbinHist(refLength, subreadLength, x_label, y_label, title, outfile, format, quantile=0.99)

def makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass'):
    title = "Reference Length vs Fraction of Reference in Alignment ({0} Subreads only)".format(label)
    x_label = "Reference Length"
    y_label = "% of Reference Aligned"

    fullPass = alnRatios[alnRatios[alnKey]]
    refLength = fullPass['RefLength']

    #readAlnLength = fullPass['rEnd'] - fullPass['rStart']
    refAlnLength = fullPass['tEnd'] - fullPass['tStart']

    alnRefRatio = refAlnLength.astype(float) / refLength.astype(float) * 100.0
    
    _makeHexbinHist(refLength, alnRefRatio, x_label, y_label, title, outfile, format, quantile=0.99)
                            

def makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refLengthRange=None):
    title = "Fraction of Reference and " + label + " Subread in Alignment"
    if refLengthRange is not None:
        title += " (ref range {0}-{1})".format(refLengthRange[0], refLengthRange[1])
    x_label = "Fraction of Reference in Alignment"
    y_label = "Fraction of " + label + " Subread in Alignment"

    if refLengthRange is None:
        fullPass = alnRatios[alnRatios[alnKey]]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]
        
    readAlnLength = fullPass['rEnd'] - fullPass['rStart']
    refAlnLength = fullPass['tEnd'] - fullPass['tStart']
    refLength = fullPass['RefLength']
    insLength = fullPass['iEnd'] - fullPass['iStart']

    alnRefRatio = refAlnLength.astype(float) / refLength.astype(float)
    alnInsRatio = readAlnLength.astype(float) / insLength.astype(float)
    
    _makeHexbinHist(alnRefRatio, alnInsRatio, x_label, y_label, title, outfile, format, quantile=0.99)

def _makeHexbinHist(x, y, x_label, y_label, title, outfile, format, quantile=None):
    nullfmt = NullFormatter()
    
    left, width = 0.1, 0.6
    bottom, height = 0.1, 0.6
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(dpi=300, figsize=(8,8))
    fig.suptitle(title)

    axHexbin = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    if quantile is not None:
        mask = n.all([x < mstats.mquantiles(x, [quantile])[0], y < mstats.mquantiles(y, [quantile])[0]], axis=0)
        x = x[mask]
        y = y[mask]
        
    x_min = n.min(x)
    x_max = n.max(x)
    y_min = n.min(y)
    y_max = n.max(y)

    axHexbin.hexbin(x, y, bins='log', edgecolors='none', cmap=plt.cm.hot)
    axHexbin.set_xlim(x_min, x_max)
    axHexbin.set_xlabel(x_label)
    axHexbin.set_ylim(y_min, y_max)
    axHexbin.set_ylabel(y_label)
    

    axHistx.hist(x, bins=100)
    axHistx.set_xlim(x_min, x_max)
    axHisty.hist(y, bins=100, orientation='horizontal')
    axHisty.set_ylim(y_min, y_max)
    for label in axHisty.get_xticklabels():
        label.set_rotation('vertical')

    fig.savefig(outfile, format=format)

def _fn(dir, pref, plot_type, suf):
    return os.path.join(dir, pref + "_" + plot_type + suf)

 
def write_summary_page(pdf_filename, args, inserts, alnRatios):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph
    
    doc = SimpleDocTemplate(pdf_filename, pagesize=letter)
    # container for the 'Flowable' objects
    styles = getSampleStyleSheet()
    elements = []
    
    input = Paragraph("""
    Input directory:
    {0}
    <br/>
    <br/>
    Files:
    """.format(os.path.abspath(args.job_directory)), styles['Normal'])
    elements.append(input)
    
    inFOFN = os.path.join(args.job_directory, "input.fofn")
    for file in open(inFOFN):
        elements.append(Paragraph(os.path.basename(file), styles['Normal']))
        
    data = []
    
    total  = 0
    fullpass = 0
    longest = 0
    AT = 0
    seq_ZMWs = 0
    for ins in inserts.itervalues(): 
        total += len(ins)
        fullpass += len(ins[ins['IsFullPass']])
        longest += len(ins[ins['IsLongest']])
        AT += len(ins[ins['IsAT']])
        seq_ZMWs += len(n.unique(ins['HoleNumber']))
                     
    elements.append(Paragraph("""
    <br/>
    Number of sequencing ZMWs: {0}
    <br/>
    Subread alignment summary:
    <br/>
    """.format(seq_ZMWs), styles['Normal']))
        
    data.append(["Subread Type", "Original", "Aligned", "Unique RefIDs aligned to"])
    data.append(["Total", total, len(alnRatios), len(set(alnRatios[:]['RefID']))])
    data.append(["FullPass", fullpass, len(alnRatios[alnRatios['IsFullPass']]), len(set(alnRatios[alnRatios['IsFullPass']]['RefID']))])
    data.append(["Longest", longest, len(alnRatios[alnRatios['IsLongest']]), len(set(alnRatios[alnRatios['IsLongest']]['RefID']))])
    data.append([SeenName, AT, len(alnRatios[alnRatios['IsAT']]), len(set(alnRatios[alnRatios['IsAT']]['RefID']))])
             
    t=Table(data)
    t.setStyle(TableStyle([('BACKGROUND', (0,0), (0,-1), colors.gray),
                           ('BACKGROUND', (0,0), (-1,0), colors.gray),
                           ('INNERGRID', (0,0), (-1,-1), .25, colors.black),
                           ('BOX', (0,0), (-1,-1), .25, colors.black)]))
    elements.append(t)
    # write the document to disk
    doc.build(elements)
    
def make_MovieDict(cmpH5):
    """
    Return a dict of MovieID --> MovieName
    """
    return cmpH5['/MovieInfo'].asDict('ID', 'Name')

def make_RefDict(cmpH5):
    """
    Return a dict of RefID --> RefName
    """
    return cmpH5['/RefInfo'].asDict('ID', 'FullName')

if __name__ == "__main__":
    global SeenName
    SeenName = "5'-3'"
    parser = argparse.ArgumentParser(description='Create some plots for transcript analyses.')
    parser.add_argument('job_directory')
    parser.add_argument('-d', '--output_directory')
    parser.add_argument("-m", '--primer_match_file')
    parser.add_argument('-p', '--output_prefix')
    parser.add_argument('-t', '--output_type')
    parser.add_argument('--read_pickle')
    parser.add_argument("--ref_size", default=None)
    parser.add_argument("--refStrandPickle", default=None)
    args = parser.parse_args()
    
    ref_size = None
    if args.ref_size is not None:
        ref_size = eval(args.ref_size)

    inFOFN = os.path.join(args.job_directory, "input.fofn")
    cmph5FN = os.path.join(args.job_directory, "data", "aligned_reads.cmp.h5")
    print "Creating cmph5 object"
    cmpH5 = factory.create(cmph5FN)
    print "Calculating reference lengths"
    refLengths = getReferenceLengths(cmpH5)

    if args.read_pickle:
        stuff = cPickle.load(open(args.read_pickle, 'rb'))
        if type(stuff) is dict:
            inserts = stuff['inserts']
            alnRatios = stuff['alns']
            refDict = stuff['RefDict']
        else:
            inserts, alnRatios = stuff
    else:    
        print "Reading inserts from input.fofn"
        inserts = getInsertsFromFofn(inFOFN, args.primer_match_file)
        print "Gathering alignment lengths"
        alnRatios = getAlignedLengthRatios(cmpH5, inserts)
        refDict = make_RefDict(cmpH5)
        cPickle.dump({'inserts':inserts,'alns':alnRatios,'MovieDict':make_MovieDict(cmpH5),'RefDict':refDict}, open(os.path.join(args.output_directory, args.output_prefix + ".pkl"), 'wb'))
        
    # make refStrandDict if needed
    refStrandDict = None
    if args.refStrandPickle is not None:
        print "Making refStrandDict"
        with open(args.refStrandPickle) as f:
            transcript_info = cPickle.load(f)
        refStrandDict = {}
        for refid, refname in refDict.iteritems():
            if refname.find('|') > 0: refname = refname[:refname.find('|')]
            refStrandDict[refid] = transcript_info[refname]['strand']

    write_summary_page(os.path.join(args.output_directory, args.output_prefix + '.summary.pdf'), args, inserts, alnRatios)
    #sys.exit(-1)

    if args.output_type == "pdf":
        pp = PdfPages(os.path.join(args.output_directory, args.output_prefix + ".figures.pdf"))
        print "Creating pdf plots"
        makeSubreadRLHistogram(alnRatios, pp, "pdf", 0.99)       
        makeFractionSubreadHistogram(alnRatios, pp, "pdf")
        makeReferenceRLHistogram(alnRatios, refLengths, pp, "pdf", 0.99)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", refStrandDict=refStrandDict)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsLongest", "Longest", refLengthRange=ref_size, refStrandDict=refStrandDict)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName, refLengthRange=ref_size, refStrandDict=refStrandDict)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", per_gene=True, refLengthRange=ref_size, refStrandDict=refStrandDict)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsLongest", "Longest", per_gene=True, refLengthRange=ref_size, refStrandDict=refStrandDict)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName, per_gene=True, refLengthRange=ref_size, refStrandDict=refStrandDict)
        makeCoverage_by_RefLength(alnRatios, pp, 'pdf')
        makeCoverage_by_RefLength(alnRatios, pp, 'pdf', 'IsLongest', 'Longest')
        makeCoverage_by_RefLength(alnRatios, pp, 'pdf', 'IsAT', SeenName)
        makeStartPositionVS(alnRatios, pp, 'pdf', refStrandDict=refStrandDict)
        makeStartPositionVS(alnRatios, pp, 'pdf', 'IsLongest', 'Longest', refStrandDict=refStrandDict)
        makeStartPositionVS(alnRatios, pp, 'pdf', 'IsAT', SeenName, refStrandDict=refStrandDict)
        makeEndPositionVS(alnRatios, pp, 'pdf', refStrandDict=refStrandDict)
        makeEndPositionVS(alnRatios, pp, 'pdf', 'IsLongest', 'Longest', refStrandDict=refStrandDict)
        makeEndPositionVS(alnRatios, pp, 'pdf', 'IsAT', SeenName, refStrandDict=refStrandDict)
        #makeRefAbundance(alnRatios, pp, 'pdf') # currently OBSOLETE
        makeRefLengthVSAbundance(alnRatios, pp, 'pdf', 'IsAT', SeenName)
        makeRefLengthVSAbundance(alnRatios, pp, 'pdf', 'IsAT', SeenName, covThreshold=.1)
        #makeFractionReferenceHistogram(alnRatios, pp, "pdf")   # not using for now, I don't think it's informative
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)
        #makeSubreadLengthVsReferenceLengthHexbinHist_LongestSubreadOnly(alnRatios, pp, "pdf", "IsAT", SeenName)
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", refLengthRange=ref_size)
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest", refLengthRange=ref_size)
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName, refLengthRange=ref_size)
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)    
        pp.close()
    elif args.output_type == "png":
        print >> sys.stderr, "Liz: made too many changes to PDF version. PNG version currently NOT supported."
        sys.exit(-1)
        makeSubreadRLHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "subread_hist", ".png"), "png", 0.99)
        makeFractionSubreadHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "frac_subread", ".png"), "png")
        makeReferenceRLHistogram(alnRatios, refLengths, _fn(args.output_directory, args.output_prefix, "ref_hist", ".png"), "png", 0.99)
        makeAlignmentPercentileDistribution(alnRatios, _fn(args.output_directory, args.output_prefix, "FPaln_start_end", ".png"), "png")
        makeAlignmentPercentileDistribution(alnRatios, _fn(args.output_directory, args.output_prefix, "LGaln_start_end", ".png"), "png", "IsLongest", "Longest")
        makeFractionReferenceHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "frac_ref", ".png"), "png")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "FPsubread_vs_ref", ".png"), "png")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "LGsubread_vs_ref", ".png"), "png", "IsLongest", "Longest")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "FPfrac_subread_vs_ref", ".png"), "png")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "LGfrac_subread_vs_ref", ".png"), "png", "IsLongest", "Longest")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "FPreflen_vs_fracref", ".png"), "png")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "LGreflen_vs_fracref", ".png"), "png", "IsLongest", "Longest")
                                                
    # concatenate the summary & figure pdf!!!
    os.system("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile={0}/{1}.pdf {0}/{1}.summary.pdf {0}/{1}.figures.pdf".format(args.output_directory, args.output_prefix))
