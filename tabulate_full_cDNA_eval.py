import os, sys, re, fnmatch, subprocess
from csv import DictReader
from pbcore.io import BasH5Reader
from collections import defaultdict
from Bio import SeqIO

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
    p1s, zmws = extract_P1(run_dir)
    p1s.sort()
    return {'machines': "|".join(machines), "cells":str(count_cells), "loading":"|".join(p1s), "seqZmws":str(sum(int(x) for x in zmws))}

def extract_P1(run_dir):
    """
    See if <run_dir>/loading_P1.txt exists, if not, produce it
    Return a list of P1 (per-movie)
    """
    file = os.path.join(run_dir, "loading_P1.txt")
    if not os.path.exists(file):
        cmd = "calc_loading_P1.py " + run_dir
        assert subprocess.check_call(cmd, shell=True) == 0
    p1s = []
    zmws = []
    with open(file, 'r') as f:
        for r in DictReader(f, delimiter=','):
            p1s.append(r['P1'])
            zmws.append(r['seqZmw'])
            
    return p1s, zmws

def get_avg_fl_len(filename):
    tmp = [len(r.seq) for  r in SeqIO.parse(open(filename),'fasta')]
    return sum(tmp)*1./len(tmp)

def read_evaled_summary(filename):
    """
    NOTE: this is now updated for ReadsOfInsert

    Number of RoI reads: 111625
    Avg. RoI readlength: 1693
    % of artificial chimeras in RoI: 132/52811 (0.2%)
    ----- GMAP ------
    Total number of ZMWs: 52679
    Total number of unmapped: 96 (0.2%)
    Avg. coverage: 98.2%
    Avg. observed accuracy: 96.1%
    Chimera rate: 1.0%
    Chimera from overpriming: 9.9%
    ----- BLASR ------
    Number of aligned refs: 10627/64877 (16.4%)
    Number of aligned ZMWs: 51592/52679 (97.9%)
    Number of well-aligned refs: 2465/64877 (3.8%)
    Number of well-aligned ZMWs: 25374/52679 (48.2%)
    """
    #print >> sys.stderr, "processing", filename
    d = {}
    with open(filename) as f:
        while True:
            if f.tell() == os.stat(filename).st_size: break
            line = f.readline().strip()
            if line.startswith('Number of RoI reads:'):
                d['roi_yield'] = line.strip().split(': ')[1]
            elif line.startswith('Avg. RoI readlength:'):
                d['roi_avg_len'] = line.strip().split(': ')[1]
            elif line.startswith('% of artificial chimeras'):
                a, b = line.strip().split(': ')[1].split()
                a, c = a.split('/')
                d['art_chimera'] = "{0} {1}".format(a, b)
                d['fl_zmws'] = c
                d['fl_good'] = str(int(c) - int(a))
            elif line.startswith('----- GMAP ------'):
                assert f.readline().startswith('Total number of ZMWs:')
                d['gmap_roi_unmapped'] = f.readline().strip().split(': ')[1]
                d['gmap_roi_cov'] = f.readline().strip().split(': ')[1]
                d['gmap_roi_acc'] = f.readline().strip().split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera rate:')
                d['gmap_roi_chimera'] = line.split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera from overpriming: ')
                d['gmap_roi_chimera_from_ovrpriming'] = line.split(': ')[1]
            elif line.startswith('----- BLASR ------'):
                line = f.readline()
                assert line.startswith('Number of aligned refs:')
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_refs'] = "{0} {1}".format(a, b)
                line = f.readline()
                assert line.startswith('Number of aligned ZMWs:')
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_ZMWs'] = "{0} {1}".format(a, b)
                line = f.readline()
                assert line.startswith('Number of well-aligned refs:')
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_refs_well'] = "{0} {1}".format(a, b)
                line = f.readline()
                assert line.startswith('Number of well-aligned ZMWs:')
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_ZMWs_well'] = "{0} {1}".format(a, b)
                break # END OF FILE!
    return d

def read_primer_summary(filename):
    """
    Return dict with '5seen', '3seen', and '53Aseen'
    """
    d = defaultdict(lambda: -1)
    with open(filename) as f:
        assert f.readline().startswith("------ 5' primer seen")
        d['5seen'] = f.readline().strip().split()[-1][1:-1]
        assert f.readline().startswith("------ 3' primer seen")
        d['3seen'] = f.readline().strip().split()[-1][1:-1]
        assert f.readline().startswith("------ 5'&3' primer seen")
        f.readline()
        assert f.readline().startswith("------ 5'&3'&polyA primer seen")
        d['53Aseen'] = f.readline().strip().split()[-1][1:-1]
    return d


def collect_info_for_name(name):
    """
    Main collecting function for runs/<name>
    Skip if primer_match/ does not exist
    """
    evaled_txt = os.path.join('primer_match', name, 'evaled_summary.txt')
    if not os.path.exists(evaled_txt):
        print >> sys.stderr, "{0} does not exist. SKIP!".format(evaled_txt)
        return None
    d = defaultdict(lambda: 'NA')
    d['name'] = name
    d['roi_primer'] = read_primer_summary(os.path.join('primer_match', name, 'reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.summary'))
    d.update(read_evaled_summary(evaled_txt))
    d['fl_avg_len'] = str(get_avg_fl_len(os.path.join('primer_match', name, 'reads_of_insert.53Aseen_trimmed_changeid.fa.non_chimera.fa')))
    d.update(extract_run_info(os.path.join('runs', name)))
    return d

#fields = ['machines', 'cells', 'ZMWs', 'P1', 'roi_yield', 'roi_avg_len', 'nonroi_avg_len', '
def write_header(f):
    f.write("Sample,Machine,Cells,SequencingZMWs,Productivity (%),RoI yield,Avg. RoI readlength,,")
    f.write("5' primer seen,3' primer seen,FL %,")
    f.write("# of FL ZMWs,Avg. FL len,Artificial Chimeras,,")
    f.write("Input # of ZMWs,# unmapped,Avg. GMAP coverage,Avg. GMAP accuracy,GMAP chimera rate,GMAP chimera from overpriming,,")
    f.write("# of aligned refs,# of aligned ZMWs,# of well-aligned refs,# of well-aligned ZMWs\n")

def write_line(f, d):
    f.write(d['name'] + ',')
    f.write(d['machines'] + ',')
    f.write(d['cells'] + ',')
    f.write(d['seqZmws'] + ',')
    f.write(d['loading'] + ',')
    f.write(d['roi_yield'] + ',')
    f.write(d['roi_avg_len'] + ',')

    f.write(',')
    f.write(d['roi_primer']['5seen'] + ',')
    f.write(d['roi_primer']['3seen'] + ',')
    f.write(d['roi_primer']['53Aseen'] + ',')
    f.write(d['fl_zmws'] + ',')
    f.write(d['fl_avg_len'] + ',')
    f.write(d['art_chimera'] + ',')
    f.write(',')

    f.write(d['fl_good'] + ',')
    f.write(d['gmap_roi_unmapped'] + ',')
    f.write(d['gmap_roi_cov'] + ',')
    f.write(d['gmap_roi_acc'] + ',')
    f.write(d['gmap_roi_chimera'] + ',')
    f.write(d['gmap_roi_chimera_from_ovrpriming'] + ',')

    f.write(',')
    f.write(d['blasr_num_refs'] + ',')
    f.write(d['blasr_num_ZMWs'] + ',')
    f.write(d['blasr_num_refs_well'] + ',')
    f.write(d['blasr_num_ZMWs_well'] + '\n')

def main():
    import fnmatch
    f = open('cDNA_summary.table.txt','w')
    write_header(f)
    for x in fnmatch.filter(os.listdir('runs'), '*'):
        name = os.path.join('runs', x)
        if not os.path.isdir(name): continue
        print >> sys.stderr, "processing", name
        d = collect_info_for_name(x)
        if d is not None:
            write_line(f, d)
    f.close()

if __name__ == "__main__":
    main()
