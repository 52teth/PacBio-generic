import os, sys, re, fnmatch, subprocess
from csv import DictReader
from pbcore.io import BasH5Reader
from collections import defaultdict
from Bio import SeqIO
from summarize_primer_info import summarize_primer_info

def extract_run_info(run_dir):
    """
    (1) Extract the machine info from the .metadata.xml file
    ex: m130524_195336_42142_c100502960310000001823071108081310_s1_p0.metadata.xml
    comes from (42)142

    (2) Extract the number of cells
    """
    machines = set()
    count_cells = 0
    for d in os.listdir(run_dir):
        d1 = os.path.join(run_dir, d)
        if not os.path.isdir(d1): continue
        for file in fnmatch.filter(os.listdir(d1), '*.metadata.xml'):
            count_cells += 1
            machine = file.split('_')[2]
            try:
                int(machine) # if possible, extract last three digit only
                machine = machine[-3:]
            except:
                pass # use the whole, probably is 'sherri'
            machines.add(machine)
    # get P1
    p1s, zmws, avgLens = extract_P1(run_dir)
    p1s.sort()
    return {'machines': "|".join(machines), "cells":str(count_cells), "loading":"|".join(map("{0:.1f}".format,p1s)), "seqZmws":str(sum(int(x) for x in zmws)), "avgLens":"|".join(map(str,avgLens))}

def extract_P1(run_dir):
    """
    See if <run_dir>/loading_P1.txt exists, if not, produce it
    Return a list of P1 (per-movie)
    """
    file = os.path.join(run_dir, "loading_P1.txt")
    if not os.path.exists(file):
        cmd = "calc_loading_P1.py " + run_dir
        assert subprocess.check_call(cmd, shell=True) == 0
    p1s = defaultdict(lambda: [])
    zmws = defaultdict(lambda: [])
    avgLens = defaultdict(lambda: [])
    names = set()
    with open(file, 'r') as f:
        for r in DictReader(f, delimiter=','):
            names.add(r['Name'])
            p1s[r['Name']].append(float(r['P1']))
            zmws[r['Name']].append(int(r['seqZmw']))
            avgLens[r['Name']].append(int(r['avgLen'])) 

    p1s = [sum(p1s[n])*1./len(p1s[n]) for n in names ]
    zmws = [sum(zmws[n]) for  n in names ]
    # the proper avgLen should be weighted by # of ZMWs in that split bax.h5 file but I'm too lazy ^_^///
    avgLens = [sum(avgLens[n])*1/len(avgLens[n]) for n in names ]
    return p1s, zmws, avgLens


def collect_info_for_name(name):
    """
    Main collecting function for runs/<name>
    Skip if eval/ does not exist
    """
    d = defaultdict(lambda: 'NA')
    d['name'] = name
    d.update(extract_run_info(os.path.join('runs', name)))
    if not os.path.exists(os.path.join('smrtpipe', name, 'isoseq_draft.primer_info.csv')):
        print >> sys.stderr, "skipping", name
        return d
    d.update(summarize_primer_info(os.path.join('smrtpipe', name, 'isoseq_draft.primer_info.csv')))
    return d

#fields = ['machines', 'cells', 'ZMWs', 'P1', 'roi_yield', 'roi_avg_len', 'nonroi_avg_len', '
def write_header(f):
    f.write("Sample,Machine,Cells,SequencingZMWs,Productivity (%),Avg. raw readlength,")
    f.write("# of CCS,5' primer seen,3' primer seen,FL %,")
    f.write("# of FL ZMWs,Avg. FL len,FL length range,Artificial Chimeras\n")

def write_line(f, d):
    f.write(d['name'] + ',')
    f.write(d['machines'] + ',')
    f.write(d['cells'] + ',')
    f.write(d['seqZmws'] + ',')
    f.write(d['loading'] + ',')
    f.write(d['avgLens'] + ',')

    f.write(str(d['num_ccs']) + ',')
    f.write(str(d['5seen']) + ',')
    f.write(str(d['3seen']) + ',')
    f.write(str(d['fl_percent']) + ',')
    f.write(str(d['fl_zmws']) + ',')
    f.write(str(d['fl_avg_len']) + ',')
    f.write(str(d['flnc_len_range']) + ',')
    f.write(d['art_chimera'])
    f.write('\n')


def main():
    import fnmatch
    pm_dict = {}
    f = open('cDNA_summary.table.txt','w')
    write_header(f)
    for x in fnmatch.filter(os.listdir('runs'), '*'):
        name = os.path.join('runs', x)
        if not os.path.isdir(name): continue
        print >> sys.stderr, "processing", name
        d = collect_info_for_name(x)
        if d is not None:
            write_line(f, d)
            pm_dict[name] = d['primer_counts']
    f.close()

    # write primer information
    with open('cDNA_summary.primer_counts.txt', 'w') as f:
        keys = set()
        for v in pm_dict.itervalues():
            if type(v) is dict:
                keys = keys.union(v.keys())
        keys = list(keys)
        keys.sort()

        f.write("sample,total," + ",".join(keys) + '\n')
        for sample, pm_count in pm_dict.iteritems():
            if type(pm_count) is not dict: continue
            f.write(os.path.basename(sample) + ',')
            total = sum(pm_count.itervalues())
            f.write(str(total))
            for k in keys:
                if k not in pm_count: f.write(",0%")
                else: f.write(",{0:.0f}%".format(pm_count[k]*100./total))
            f.write('\n')


if __name__ == "__main__":
    main()
