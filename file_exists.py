#!/usr/bin/env python
import os, sys

file = sys.argv[1]
for line in open(file):
    if not os.path.exists(line.strip()):
        print line.strip(), "does not exist"
