#!/usr/bin/env python
"""
Quick script to go through each subfolder in runs/
ex: H1_1to2k/, H1_2to3k/
and display a list of SMRT cell samples associated with it.

ex:
=========
H1_1to2k
=========
<path> - <WellName> - <CollectionNumber> - <Name>
"""

import os, sys, fnmatch
from xml.dom.minidom import parse

def get_cell_info(xml_filename):
    dom = parse(xml_filename)
    wellname = dom.getElementsByTagName('WellName')[0].firstChild.nodeValue
    collection = dom.getElementsByTagName('CollectionNumber')[0].firstChild.nodeValue
    samplename = dom.getElementsByTagName('Sample')[0].getElementsByTagName('Name')[0].firstChild.nodeValue
    comments = dom.getElementsByTagName('Comments')[0].firstChild.nodeValue
    return wellname, collection, samplename, comments

def main():
    for dirname in os.listdir('runs/'):
        if not os.path.isdir('runs/'+dirname): continue
        print "========================================="
        print dirname
        print "========================================="
        rnames = os.listdir(os.path.join('runs/', dirname))
        rnames.sort()
        for runname in rnames:
            xml =  os.path.join('runs/', dirname, runname, fnmatch.filter(os.listdir(os.path.join('runs/', dirname, runname)), '*.metadata.xml')[0])
            print runname+":",os.path.basename(xml)[:-19],"-".join(get_cell_info(xml))
            
        
if __name__ == "__main__":
    main()
        
