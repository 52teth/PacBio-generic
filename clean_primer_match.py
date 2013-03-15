#!/usr/bin/env python
import os, sys, fnmatch, shutil


for d0 in os.listdir('.'):
    if not os.path.isdir(d0): continue
    for d1 in os.listdir(d0):
        if d1.startswith('primer_match'):
            for d2 in os.listdir(os.path.join(d0,d1)):
                x2 = os.path.join(d0,d1,d2)
                if not os.path.isdir(x2): continue
                for d3 in fnmatch.filter(os.listdir(x2), 'output*'):
                    x3 = os.path.join(x2, d3)
                    if not os.path.isdir(x3): continue
                    for file in os.listdir(x3):
                        if file!='all.RESULT.reduceMax':
                            x4 = os.path.join(x3, file)
                            print >> sys.stderr, "removing", x4
                            if os.path.isdir(x4):
                                shutil.rmtree(x4)
                            else:
                                os.remove(x4)

                
