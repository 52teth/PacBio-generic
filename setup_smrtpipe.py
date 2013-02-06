#!/usr/bin/env python
"""
Sets up directories in smrtpipe/ according to runs/ directory
"""
import os, sys, subprocess

for name in os.listdir('runs/'):
    d_run = os.path.join('runs', name)
    if not os.path.isdir(d_run): continue

    d = os.path.join('smrtpipe/', name)
    if not os.path.exists(d):
        print >> sys.stderr, "creating directory", d
        os.makedirs(d)

    cmd = "find {0}/*/Analysis_Results/*.bas.h5 > {1}/input.fofn".format(os.path.abspath(d_run),d)
    if os.system(cmd)!=0:
        print >> sys.stderr, "problem running", cmd
        continue

        

