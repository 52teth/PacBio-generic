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
    p1s, zmws, avgLens = extract_P1(run_dir)
    p1s.sort()
    return {'machines': "|".join(machines), "cells":str(count_cells), "loading":"|".join(p1s), "seqZmws":str(sum(int(x) for x in zmws)), "avgLens":"|".join(avgLens)}

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
    avgLens = []
    with open(file, 'r') as f:
        for r in DictReader(f, delimiter=','):
            p1s.append(r['P1'])
            zmws.append(r['seqZmw'])
            avgLens.append(r['avgLen']) 
    return p1s, zmws, avgLens

def read_evaled_summary(filename):
    """
    NOTE: this is now updated for ReadsOfInsert

    RoI length distribution: 1000-2000
    Avg. RoI length: 1290
    flnc length distribution: 1000-2000
    Estimated artificial concatemer rate: 0/38502 (0.00%)
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
            if line.startswith('RoI length distribution: '):
                d['roi_len_range'] = line.strip().split(': ')[1]
            elif line.startswith('Avg. RoI length: '):
                d['avg_roi_len'] = int(line.strip().split(': ')[1])
            elif line.startswith('flnc length distribution: '):
                d['flnc_len_range'] = line.strip().split(': ')[1]
            elif line.startswith('Estimated artificial concatemer rate: '):
                a, b = line.strip().split(': ')[1].split()
                a, c = a.split('/')
                d['art_chimera'] = "{0} {1}".format(a, b)
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
                d['gmap_roi_chimera_from_overpriming'] = line.split(': ')[1]
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

def read_classify_summary(filename):
    """
    Number of reads of insert=155623
    Number of five prime reads=108137
    Number of three prime reads=112952
    Number of poly-A reads=112937
    Number of filtered short reads=6320
    Number of non-full-length reads=56789
    Number of full-length reads=92514
    Number of full-length non-chimeric reads=92463
    Average full-length non-chimeric read length=1718
    """
    d = defaultdict(lambda: -1)
    with open(filename) as f:
        line = f.readline().strip()
        assert line.startswith("Number of reads of insert=")
        d['roi_yield'] = int(line.split('=')[1])
        line = f.readline().strip()
        assert line.startswith("Number of five prime reads=")
        d['5seen'] = int(line.split('=')[1]) * 100. / d['roi_yield']
        line = f.readline().strip()
        assert line.startswith("Number of three prime reads=")
        d['3seen'] = int(line.split('=')[1]) * 100. / d['roi_yield']
        line = f.readline().strip()
        assert line.startswith("Number of poly-A reads=")
        line = f.readline().strip()
        assert line.startswith('Number of filtered short reads')
        line = f.readline().strip()
        assert line.startswith('Number of non-full-length reads')
        line = f.readline().strip()
        assert line.startswith('Number of full-length reads')
        num_fl = int(line.split('=')[1])
        line = f.readline().strip()
        assert line.startswith("Number of full-length non-chimeric reads=")
        d['fl_zmws'] = int(line.split('=')[1])
        num_chim = num_fl - d['fl_zmws']
        d['art_chimera'] = "{0} ({1:.2f}%)".format(num_chim, num_chim*100./num_fl)
        d['53Aseen'] = d['fl_zmws'] * 100. / d['roi_yield']
        line = f.readline().strip()
        assert line.startswith("Average full-length non-chimeric read length")
        d['fl_avg_len'] = int(line.split('=')[1])
    return d


def collect_info_for_name(name):
    """
    Main collecting function for runs/<name>
    Skip if eval/ does not exist
    """
    d = defaultdict(lambda: 'NA')
    d['name'] = name
    d.update(read_classify_summary(os.path.join('smrtpipe', name, 'data', 'classify_summary.txt')))
    evaled_txt = os.path.join('eval', name, 'evaled_summary.txt')
    if os.path.exists(evaled_txt):
        d.update(read_evaled_summary(evaled_txt))
    else:
        print >> sys.stderr, "{0} does not exist. SKIP!".format(evaled_txt)
        roi_fasta = os.path.join('smrtpipe', name, 'data', 'reads_of_insert.fasta')
        roi_lens = [len(r.seq) for r in SeqIO.parse(open(roi_fasta), 'fasta')]
        d['avg_roi_len'] = sum(roi_lens)*1./len(roi_lens)
    d.update(extract_run_info(os.path.join('runs', name)))
    return d

#fields = ['machines', 'cells', 'ZMWs', 'P1', 'roi_yield', 'roi_avg_len', 'nonroi_avg_len', '
def write_header(f):
    f.write("Sample,Machine,Cells,SequencingZMWs,Productivity (%),Avg. raw readlength,RoI yield,Avg. RoI readlength,RoI readlength range,,")
    f.write("5' primer seen,3' primer seen,FL %,")
    f.write("# of FL ZMWs,Avg. FL len,FL length range,Artificial Chimeras,,")
    f.write("Input # of ZMWs,# unmapped,Avg. GMAP coverage,Avg. GMAP accuracy,GMAP chimera rate,GMAP chimera from overpriming,,")
    f.write("# of aligned refs,# of aligned ZMWs,# of well-aligned refs,# of well-aligned ZMWs\n")

def write_line(f, d):
    f.write(d['name'] + ',')
    f.write(d['machines'] + ',')
    f.write(d['cells'] + ',')
    f.write(d['seqZmws'] + ',')
    f.write(d['loading'] + ',')
    f.write(d['avgLens'] + ',')
    f.write(str(d['roi_yield']) + ',')
    f.write(str(d['avg_roi_len']) + ',')
    f.write(str(d['roi_len_range']) + ',')

    f.write(',')
    f.write(str(d['5seen']) + ',')
    f.write(str(d['3seen']) + ',')
    f.write(str(d['53Aseen']) + ',')
    f.write(str(d['fl_zmws']) + ',')
    f.write(str(d['fl_avg_len']) + ',')
    f.write(str(d['flnc_len_range']) + ',')
    f.write(d['art_chimera'] + ',')
    f.write(',')

    f.write(str(d['fl_zmws']) + ',')
    f.write(d['gmap_roi_unmapped'] + ',')
    f.write(d['gmap_roi_cov'] + ',')
    f.write(d['gmap_roi_acc'] + ',')
    f.write(d['gmap_roi_chimera'] + ',')
    f.write(d['gmap_roi_chimera_from_overpriming'] + ',')

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
