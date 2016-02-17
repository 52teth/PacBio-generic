#!/usr/bin/env python

import numpy as np
from csv import DictReader
from collections import defaultdict


"""
ex:
id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera
m151104_115017_42133_c100929862550000001823209905251602_s1_p0/18/1991_72_CCS,-,1,1,1,39,1958,1989,3,0
m151104_115017_42133_c100929862550000001823209905251602_s1_p0/35/1576_73_CCS,-,1,1,1,32,1535,1567,2,0
m151104_115017_42133_c100929862550000001823209905251602_s1_p0/43/31_1527_CCS,+,1,1,1,31,1527,1560,0,0
m151104_115017_42133_c100929862550000001823209905251602_s1_p0/48/31_1656_CCS,+,1,1,1,31,1656,1685,3,0
"""

def summarize_primer_info(filename="isoseq_draft.primer_info.csv"):
    num_read, num_read_fiveseen, num_read_threeseen, num_read_5threeseen, num_read_53Aseen, num_read_FLNC = 0, 0, 0, 0, 0, 0
    num_chim = 0
    pm_count = defaultdict(lambda: 0) # primer --> count
    flnc_lengths = []
    with open(filename, 'r') as f:
        for r in DictReader(f, delimiter=','):
            see5 = r['fiveseen']=='1'
            see3 = r['threeseen']=='1'
            seeA = r['polyAseen']=='1'
            nonchim = r['chimera']=='0'
            if r['primer']!='NA':
                pm_count[r['primer']] += 1

            num_read += 1
            num_read_fiveseen += see5
            num_read_threeseen += see3
            num_read_5threeseen += see5 and see3
            num_read_53Aseen += see5 and see3 and seeA
            num_read_FLNC += see5 and see3 and seeA and nonchim
            if see5 and see3 and seeA:
                num_chim += not nonchim
            if see5 and see3 and seeA and nonchim:
                flnc_lengths.append(int(r['polyAend'])-int(r['fiveend']))

        flnc_lengths = np.array(flnc_lengths)
        _low, _high = np.percentile(flnc_lengths, [5, 95])

    return {'num_ccs': num_read, '5seen':num_read_fiveseen, '3seen':num_read_threeseen, \
            '53seen': num_read_5threeseen, '53Aseen': num_read_53Aseen, 'fl_zmws': num_read_FLNC, \
            'primer_counts': dict(pm_count), 'flnc_avg_len': np.mean(flnc_lengths), \
            'flnc_len_range': "{0:.0f}-{1:.0f}".format(_low, _high), \
            "art_chimera": "{0} ({1:.2f}%)".format(num_chim, num_chim*100./num_read_FLNC)}

if __name__ == "__main__":

    d = summarize_primer_info()
    print "# of reads:", d['num_ccs']
    print "# of 5' reads:", d['5seen']
    print "# of 3' reads:", d['3seen']
    print "# of 5'&3' reads:", d['53seen']
    print "# of 5'&3'&polyA reads:", d['53Aseen']
    print "# of FLNC reads:", d['fl_zmws']
    pm_count = d['primer_counts']
    print "------ Primer Match breakdown ----"
    keys = pm_count.keys()
    keys.sort()
    total = sum(pm_count.itervalues())
    for k in keys:
        print "F{0}/R{0}: {1} ({2:.1f}%)".format(k, pm_count[k], pm_count[k]*100./total)
