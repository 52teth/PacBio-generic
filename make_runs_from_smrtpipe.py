#!/usr/bin/env python

"""
Quick script that takes the symbolic links from smrtpipe/
And make up the corressponding runs/
"""

import os,glob
from collections import defaultdict

global_seen = set()

dirs = os.listdir('smrtpipe')
for d in dirs:
    fofn = os.path.join('smrtpipe', d, 'input.fofn')
    seen = defaultdict(lambda: [])
    for line in open(fofn):
        line =line.strip()
        if not line.endswith('.1.bax.h5'): continue
        if line in global_seen:
            raise Exception, "{0} is seen twice! NOT OK!".format(line)
        global_seen.add(line)
        x = os.path.dirname(line)
        assert x.endswith('Analysis_Results')
        path = os.path.dirname(x)
        name = path.split('/')[-1]
        seen[name].append(path)
    d2 = os.path.join('runs', d)
    if not os.path.exists(d2): os.makedirs(d2)
    os.chdir(d2)
    for name,links in seen.iteritems():
        for i,link in enumerate(links):
            os.system("ln -s {0} {1}_{2}".format(link, name, i))
    os.chdir('../../')
    print d, "done"
    
