#!/usr/bin/env python

import os
import numpy as n
from re import sub
import argparse
import pickle

from pbcore.io.BasH5IO import BasH5
from pbcore.io.cmph5 import factory, alignmentPairMap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullFormatter, LinearLocator

from scipy.stats import mode, mstats

def getInsertsFromBasH5(bash5FN):
    bash5 = BasH5(bash5FN)
    rgnTable = bash5.rgnTable
    data = []

    good, bad = 0,0
    dummy = 1
    if 'Adapter' in bash5.rgnTable.rgnDS.attrs.get('RegionTypes',[]):
        for hn in bash5.getSequencingZMWs():
            print dummy
            dummy += 1
            # get indices first
            adapters = rgnTable.getAdapterRegionForZMW(hn)
            hqregion = rgnTable.getHQRegionForZMW(hn)
            inserts = rgnTable.getInsertRegionForZMW(hn)

            adapterStarts = [n.min(x) for x in adapters]
            adapterEnds = [n.max(x) for x in adapters]
            hqStart = hqregion[0][0]
            hqEnd = hqregion[0][1]

            for i in inserts:
                #print "Examining", i
                iStart, iEnd = list(i)


                is_in_hq = hqStart <= iStart < hqEnd or hqStart <= iEnd < hqEnd
                if not is_in_hq:
                    continue

#                # check if any part of insert in HQRegion
#                hqRange = range(hqStart, hqEnd)
#                if not iStart in hqRange and not iEnd in hqRange:
#                    #print "Insert not in HQRange", hqregion
#                    continue

                 # check if insert is full pass
                isFullPass = iStart in adapterEnds and iEnd in adapterStarts

                isHQTrimmed = False
                # trim by HQRegion
                if iStart < hqStart:
                    iStart = hqStart
                    isHQTrimmed = True
                if iEnd > hqEnd:
                    iEnd = hqEnd
                    isHQTrimmed = True

                insert = (hn, iStart, iEnd, isFullPass, isHQTrimmed)
                #print "Adding", insert
                data.append(insert)
    return n.array(data, dtype=[('HoleNumber', '<i4'), ('rStart', '<i4'), ('rEnd', '<i4'), ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool)])

def getInsertsFromFofn(inputFOFN):
    files = [x.strip() for x in open(inputFOFN, 'r')]
    inserts = {}
    for f in files:
        print "Reading", f
        movieName = sub('.pls.h5|.bas.h5', '', os.path.basename(f))
        inserts[movieName] = getInsertsFromBasH5(f)

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

                for a in alns:
                    rStart = a['rStart']
                    rEnd = a['rEnd']
                    tStart = a['tStart']
                    tEnd = a['tEnd']
                    refGroupID = a['RefGroupID']
                    refLength = refLengthDict[refIdDict[refGroupID]]

                    for i in inss:
                        iStart = i['rStart']
                        iEnd = i['rEnd']
                        iIsFullPass = i['IsFullPass']
                        iIsHQTrimmed = i['IsHQTrimmed']

                        if rStart >= iStart and rEnd <= iEnd:
                            results.append((movieID, refIdDict[refGroupID], hn, rStart, rEnd, tStart, tEnd, iStart, iEnd, refLength, iIsFullPass, iIsHQTrimmed))

    return n.array(results, dtype=[('MovieID', '<i2'), ('RefID', '<i4'), ('HoleNumber', '<i4'),
                                   ('rStart', '<i4'), ('rEnd', '<i4'),
                                   ('tStart', '<i4'), ('tEnd', '<i4'),
                                   ('iStart', '<i4'), ('iEnd', '<i4'), ('RefLength', '<i4'),
                                   ('IsFullPass', n.bool), ('IsHQTrimmed', n.bool)])

def makeFractionSubreadHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Subread in Alignment Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]

    for l, label in zip((alnRatios, HQRegion, fullPass), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads")):
        alnLength = l['rEnd'] - l['rStart']
        srLength = l['iEnd'] - l['iStart']
        alnSubRatio = alnLength.astype(float) / srLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper center')
    fig.savefig(outfile, format=format)

def makeSubreadRLHistogram(alnRatios, outfile, format, quantile=None):
    fig = plt.figure(dpi=300, figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Aligned Subread RL Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]

    for l, label in zip((alnRatios, HQRegion, fullPass), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads")):
        srLength = l['iEnd'] - l['iStart']

        if not quantile == None:
            srLength = srLength[srLength < mstats.mquantiles(srLength, [quantile])[0]]

        ax.hist(srLength, bins=100, histtype='step', label=label)

    ax.legend(loc='upper center')
    fig.savefig(outfile, format=format)

def makeFractionReferenceHistogram(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Fraction of Reference in Alignment Histogram")

    fullPass = alnRatios[alnRatios['IsFullPass']]
    HQRegion = alnRatios[n.any([alnRatios['IsFullPass'], alnRatios['IsHQTrimmed']], axis=0)]
    colors = ('b', 'r', 'g')

    for l, label in zip((alnRatios, HQRegion, fullPass), ("All Subreads", "HQRegion Subreads", "HQRegion Full-Pass Subreads")):
        alnLength = l['tEnd'] - l['tStart']
        refLength = l['RefLength']
        alnSubRatio = alnLength.astype(float) / refLength.astype(float)
        ax.hist(alnSubRatio, bins=100, histtype='step', label=label)

    ax.legend(loc='upper center')
    fig.savefig(outfile, format=format)

def makeReferenceRLHistogram(alnRatios, refLengths, outfile, format, quantile=None):
    fig = plt.figure(dpi=300, figsize=(6, 6))
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
    fig.savefig(outfile, format=format)

def makeAlignmentPercentileDistribution(alnRatios, outfile, format):
    fig = plt.figure(dpi=300, figsize=(8, 8))
    ax1 = fig.add_subplot(211)
    #ax1.set_title("Normalized Reference Start/End Positions")
    ax1.set_title("Alignment Reference Start/End Positions")

    ax2 = fig.add_subplot(212)
    #ax2.set_title("Normalized Reference Coverage")
    ax2.set_title("Reference Coverage")

    fullPass = alnRatios[alnRatios['IsFullPass']]
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
                                 
    ax1.set_ylim(0, max_y * 1.1)
    ax1.set_xlim(-0.05, 1.05)
    ax1.legend(loc='upper center', prop={'size': 'small'})
    
    fig.savefig(outfile, format=format)

def makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, outfile, format):
    title = "Full-Pass Subread Length vs Reference Length"
    x_label = "Reference Length"
    y_label = "Subread Length"

    fullPass = alnRatios[alnRatios['IsFullPass']]
    subreadLength = fullPass['iEnd'] - fullPass['iStart']
    refLength = fullPass['RefLength']
    
    _makeHexbinHist(refLength, subreadLength, x_label, y_label, title, outfile, format, quantile=0.99)

def makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, outfile, format):
    title = "Reference Length vs Fraction of Reference in Alignment"
    x_label = "Reference Length"
    y_label = "% of Reference Aligned"

    fullPass = alnRatios[alnRatios['IsFullPass']]
    refLength = fullPass['RefLength']
    readAlnLength = fullPass['rEnd'] - fullPass['rStart']
    refAlnLength = fullPass['tEnd'] - fullPass['tStart']

    alnRefRatio = refAlnLength.astype(float) / refLength.astype(float) * 100.0
    
    _makeHexbinHist(refLength, alnRefRatio, x_label, y_label, title, outfile, format, quantile=0.99)
                            

def makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, outfile, format):
    title = "Fraction of Reference and Subread in Alignment"
    x_label = "Fraction of Reference in Alignment"
    y_label = "Fraction of Insert in Alignment"

    fullPass = alnRatios[alnRatios['IsFullPass']]
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create some plots for transcript analyses.')
    parser.add_argument('job_directory')
    parser.add_argument('-d', '--output_directory')
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
        alnRatios = pickle.load(open(args.read_pickle, 'rb'))

    else:
        print "Reading inserts from input.fofn"
        inserts = getInsertsFromFofn(inFOFN)
        print "Gathering alignment lengths"
        alnRatios = getAlignedLengthRatios(cmpH5, inserts)
        pickle.dump(alnRatios, open(os.path.join(args.output_directory, args.output_prefix + ".pkl"), 'wb'))

    if args.output_type == "pdf":
        pp = PdfPages(os.path.join(args.output_directory, args.output_prefix + ".pdf"))

        print "Creating pdf plots"
        makeSubreadRLHistogram(alnRatios, pp, "pdf", 0.99)
        makeFractionSubreadHistogram(alnRatios, pp, "pdf")
        makeReferenceRLHistogram(alnRatios, refLengths, pp, "pdf", 0.99)
        makeAlignmentPercentileDistribution(alnRatios, pp, "pdf")
        makeFractionReferenceHistogram(alnRatios, pp, "pdf")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, pp, "pdf")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, pp, "pdf")
    
        pp.close()
    elif args.output_type == "png":
        
        makeSubreadRLHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "subread_hist", ".png"), "png", 0.99)
        makeFractionSubreadHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "frac_subread", ".png"), "png")
        makeReferenceRLHistogram(alnRatios, refLengths, _fn(args.output_directory, args.output_prefix, "ref_hist", ".png"), "png", 0.99)
        makeAlignmentPercentileDistribution(alnRatios, _fn(args.output_directory, args.output_prefix, "aln_start_end", ".png"), "png")
        makeFractionReferenceHistogram(alnRatios, _fn(args.output_directory, args.output_prefix, "frac_ref", ".png"), "png")
        makeSubreadLengthVsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "subread_vs_ref", ".png"), "png")
        makeFractionReferencevsFractionSubreadHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "frac_subread_vs_ref", ".png"), "png")
        makeFractionReferencevsReferenceLengthHexbinHist(alnRatios, _fn(args.output_directory, args.output_prefix, "reflen_vs_fracref", ".png"), "png")
                                                
    
