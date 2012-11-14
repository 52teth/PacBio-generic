#!/usr/bin/env python

import os, sys
import numpy as n
from re import sub
import argparse
import pickle
import pdb

from pbcore.io.BasH5IO import BasH5
from pbcore.io.cmph5 import factory, alignmentPairMap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter, LinearLocator

from scipy.stats import mode, mstats

from collections import defaultdict

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
    return grouped

def getPolyAT(filename):
    """
    Read polyAT report file (from seqclean)
    Movie   HoleNumber      rStart  rEnd
    
    Returns a dict of (movieName, holeNumber) --> list of (start, end) of polyAT
    """
    from csv import DictReader
    d = defaultdict(lambda: [])
    for r in DictReader(open(filename), delimiter='\t'):
        d[(r['Movie'],r['HoleNumber'])].append((int(r['rStart']), int(r['rEnd'])))
    return d
        
def hasPolyAT(d, movie, hn, start, end):
    """
    Best used with functools.partial to create a boolean function for checking
    whether an insert is in the dictionary
    """
    k = movie, hn
    if k not in d:
        return False
    for s, e in d[k]:
        if s <= start < end <= e:
            return True
    return False    

def getInsertsFromBasH5(bash5FN, func_is_polyAT):
    bash5 = BasH5(bash5FN)
    rgnTable = bash5.rgnTable
    data = []
    

    good, bad = 0,0
    dummy = 1
    if 'Adapter' in bash5.rgnTable.rgnDS.attrs.get('RegionTypes',[]):
        for hn in bash5.getSequencingZMWs():
            #print dummy
            dummy += 1
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
            #max_insert_length = max(i[1]-i[0] for i in hqinserts)
            #max_insert_length = 0

            for i in hqinserts:
                #print "Examining", i
                iStart, iEnd = list(i)


                # LIZ: no longer needed this part as I pre-filtered to hqinserts
#                is_in_hq = hqStart <= iStart < hqEnd or hqStart <= iEnd < hqEnd
#                if not is_in_hq:
#                    continue

#                # check if any part of insert in HQRegion
#                hqRange = range(hqStart, hqEnd)
#                if not iStart in hqRange and not iEnd in hqRange:
#                    #print "Insert not in HQRange", hqregion
#                    continue

                 # check if insert is full pass
                isFullPass = iStart in adapterEnds and iEnd in adapterStarts
                
                # technically this is not always right if HQtrim happens, but oh well should be OK
                #isLongest = (iEnd-iStart) == max_insert_length
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
    import functools
    
    primer_match_dict = getPolyAT(primer_match_filename)
    
    files = [x.strip() for x in open(inputFOFN, 'r')]
    inserts = {}
    for f in files:
        print "Reading", f
        movieName = sub('.pls.h5|.bas.h5', '', os.path.basename(f))
        inserts[movieName] = getInsertsFromBasH5(f, functools.partial(hasPolyAT, primer_match_dict, movieName))

    return inserts

def getReferenceLengths(cmph5):
    return cmph5['/RefInfo'].asRecArray()['Length']


def getAlignedLengthRatios(cmph5, inserts):
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

                    results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, tStart, tEnd, iStart, iEnd, refLength, iIsFullPass, iIsHQTrimmed, iIsLongest, iIsAT))
                

#                for a in alns:
#                    rStart = a['rStart']
#                    rEnd = a['rEnd']
#                    tStart = a['tStart']
#                    tEnd = a['tEnd']
#                    refGroupID = a['RefGroupID']
#                    refLength = refLengthDict[refIdDict[refGroupID]]
#
#                    for i in inss:
#                        iStart = i['rStart']
#                        iEnd = i['rEnd']
#                        iIsFullPass = i['IsFullPass']
#                        iIsHQTrimmed = i['IsHQTrimmed']
#                        iIsLongest = i['IsLongest']
#                        iIsAT = i['IsAT']
#
#                        if rStart >= iStart and rEnd <= iEnd:
#                            results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, tStart, tEnd, iStart, iEnd, refLength, iIsFullPass, iIsHQTrimmed, iIsLongest, iIsAT))
#
    return n.array(results, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),
                                   ('rStart', '<i4'), ('rEnd', '<i4'),
                                   ('tStart', '<i4'), ('tEnd', '<i4'),
                                   ('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'),
                                   ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool), ('IsLongest', n.bool),
                                   ('IsAT', n.bool)])
                                        

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
        alnSubRatio = alnLength.astype(float) / srLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper left')
    ax.set_xlabel("Alignment length/Subread Length ratio")
    ax.set_ylabel("Count")
    fig.savefig(outfile, format=format)

def makeSubreadRLHistogram(alnRatios, outfile, format, quantile=None):
    fig = plt.figure(dpi=300, figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Subread RL Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    Longest = alnRatios[alnRatios['IsLongest']]
    AT = alnRatios[alnRatios['IsAT']]
    
    for l, label in zip((alnRatios, HQRegion, fullPass, Longest, AT), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads", "HQRegion Longest Subreads", "HQRegion " + SeenName)):
        srLength = l['iEnd'] - l['iStart']

        if not quantile == None:
            srLength = srLength[srLength < mstats.mquantiles(srLength, [quantile])[0]]

        ax.hist(srLength, bins=100, histtype='step', label=label)

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

def makeAlignmentPercentileDistribution(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass', refLengthRange=None):
    fig = plt.figure(dpi=300, figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    #ax1.set_title("Normalized Reference Start/End Positions")
    ax1.set_title("Alignment Reference Start/End Positions ({0} Subreads only)".format(label))

    ax2 = fig.add_subplot(212)
    #ax2.set_title("Normalized Reference Coverage")
    if refLengthRange is None:
        ax2.set_title("Reference Coverage")
    else:
        ax2.set_title("Reference Coverage (size:{0}-{1})".format(refLengthRange[0], refLengthRange[1]))
    
    if refLengthRange is None:
        fullPass = alnRatios[alnRatios[alnKey]]
    else:
        fullPass = alnRatios[alnRatios[alnKey]&(refLengthRange[0]<=alnRatios['RefLength'])&(alnRatios['RefLength']<=refLengthRange[1])]
        
    refLength = fullPass['RefLength']
    startPercentile = fullPass['tStart'].astype(float) / refLength.astype(float)
    endPercentile = fullPass['tEnd'].astype(float) / refLength.astype(float)

    normalized_alignments = n.concatenate([n.arange(s, e, 0.01) for s, e in zip(n.round(startPercentile, decimals=2), n.around(endPercentile, decimals=2))])

    max_y = 0
    for l, label in zip([startPercentile, endPercentile], ['Alignment Start', 'Alignment End']):
        num, bins, patches = ax1.hist(l, bins=100,
                                      histtype='step', label=label, normed=True)

        if n.max(num) > max_y:
            max_y = n.max(num)

    ax2.hist(normalized_alignments, bins=100, histtype='stepfilled')
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel("Normalized Reference Sequence Length")
    ax2.set_ylabel("Count")
                                 
    ax1.set_ylim(0, max_y * 1.1)
    ax1.set_xlim(-0.05, 1.05)
    ax1.legend(loc='upper center', prop={'size': 'small'})
    
    fig.savefig(outfile, format=format)

def makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', title='Full-Pass'):
    title = title + " Subread Length vs Reference Length"
    x_label = "Reference Length"
    y_label = "Subread Length"

    fullPass = alnRatios[alnRatios[alnKey]]
    subreadLength = fullPass['iEnd'] - fullPass['iStart']
    refLength = fullPass['RefLength']
    
    _makeHexbinHist(refLength, subreadLength, x_label, y_label, title, outfile, format, quantile=0.99)
    
def makeSubreadLengthVsReferenceLengthHexbinHist_LongestSubreadOnly(alnRatios, outfile, format, alnKey='IsFullPass', title='Full-Pass'):
    title = title + " Subread Length vs Reference Length (longest subread only!)"
    x_label = "Reference Length"
    y_label = "Subread Length"

    subreads = alnRatios[alnRatios[alnKey]]
    refIDs = list(set(subreads['RefID']))
    
    subreadLength = []
    refLength = []
    # for each refID, find the one with the max length
    for refID in refIDs:
        x = alnRatios[alnRatios['RefID']==refID]
        subreadLength.append(max(x['iEnd']-x['iStart']))
        refLength.append(x[0]['RefLength'])

    _makeHexbinHist(n.array(refLength), n.array(subreadLength), x_label, y_label, title, outfile, format, quantile=0.99)

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
                            

def makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, outfile, format, alnKey='IsFullPass', label='Full-Pass'):
    title = "Fraction of Reference and " + label + " Subread in Alignment"
    x_label = "Fraction of Reference in Alignment"
    y_label = "Fraction of " + label + " Subread in Alignment"

    fullPass = alnRatios[alnRatios[alnKey]]
    readAlnLength = fullPass['rEnd'] - fullPass['rStart']
    refAlnLength = fullPass['tEnd'] - fullPass['tStart']
    refLength = fullPass['RefLength']
    insLength = fullPass['iEnd'] - fullPass['iStart']

    alnRefRatio = refAlnLength.astype(float) / refLength.astype(float)
    alnInsRatio = readAlnLength.astype(float) / insLength.astype(float)
    
    _makeHexbinHist(alnRefRatio, alnInsRatio, x_label, y_label, title, outfile, format)

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

    if not quantile==None:
        mask = n.all([x < mstats.mquantiles(x, [quantile])[0], y < mstats.mquantiles(y, [quantile])[0]], axis=0)
        x = x[mask]
        y = y[mask]

    axHexbin.hexbin(x, y, bins='log')
    axHexbin.set_xlim(0, n.max(x))
    axHexbin.set_xlabel(x_label)
    axHexbin.set_ylim(0, n.max(y))
    axHexbin.set_ylabel(y_label)

    axHistx.hist(x, bins=100)
    axHistx.set_xlim(0, n.max(x))
    axHisty.hist(y, bins=100, orientation='horizontal')
    axHisty.set_ylim(0, n.max(y))
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
        
    elements.append(Paragraph("""
    <br/>
    <br/>
    Subread alignment summary:
    <br/>
    """, styles['Normal']))
        
    data = []
    
    total  = 0
    fullpass = 0
    longest = 0
    AT = 0
    for ins in inserts.itervalues(): 
        total += len(ins)
        fullpass += len(ins[ins['IsFullPass']])
        longest += len(ins[ins['IsLongest']])
        AT += len(ins[ins['IsAT']])
        
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
    return cmpH5['/MovieInfo'].asDict('ID', 'FullName')

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
    args = parser.parse_args()

    inFOFN = os.path.join(args.job_directory, "input.fofn")
    cmph5FN = os.path.join(args.job_directory, "data", "aligned_reads.cmp.h5")

    print "Creating cmph5 object"
    cmpH5 = factory.create(cmph5FN)
    print "Calculating reference lengths"
    refLengths = getReferenceLengths(cmpH5)
   
    if args.read_pickle:
        stuff = pickle.load(open(args.read_pickle, 'rb'))
        if type(stuff) is dict:
            inserts = stuff['inserts']
            alnRatios = stuff['alns']
        else:
            inserts, alnRatios = stuff
    else:
        print "Reading inserts from input.fofn"
        inserts = getInsertsFromFofn(inFOFN, args.primer_match_file)
        print "Gathering alignment lengths"
        alnRatios = getAlignedLengthRatios(cmpH5, inserts)
        pickle.dump({'inserts':inserts,'alns':alnRatios,'MovieDict':make_MovieDict(cmpH5),'RefDict':make_RefDict(cmpH5)}, open(os.path.join(args.output_directory, args.output_prefix + ".pkl"), 'wb'))

    write_summary_page(os.path.join(args.output_directory, args.output_prefix + '.summary.pdf'), args, inserts, alnRatios)
    #sys.exit(-1)

    if args.output_type == "pdf":
        pp = PdfPages(os.path.join(args.output_directory, args.output_prefix + ".figures.pdf"))

        print "Creating pdf plots"
        makeSubreadRLHistogram(alnRatios, pp, "pdf", 0.99)
        makeFractionSubreadHistogram(alnRatios, pp, "pdf")
        makeReferenceRLHistogram(alnRatios, refLengths, pp, "pdf", 0.99)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf")
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName)
        #makeAlignmentPercentileDistribution(alnRatios, pp, "pdf", "IsAT", SeenName, (1000,2000))
        makeFractionReferenceHistogram(alnRatios, pp, "pdf")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)
        makeSubreadLengthVsReferenceLengthHexbinHist_LongestSubreadOnly(alnRatios, pp, "pdf", "IsAT", SeenName)
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsLongest", "Longest")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf", "IsAT", SeenName)
    
        pp.close()
    elif args.output_type == "png":
        
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
