#!/usr/bin/env python
import os, sys

fofn = sys.argv[1]
for line in open(fofn):
    if not os.path.exists(line.strip()):
        print "does not exists:", line.strip()
