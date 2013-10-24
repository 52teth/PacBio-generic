import os, sys, re, fnmatch
from csv import DictReader

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
        for file in fnmatch.filter(os.listdir(d1), '*.metadata.xml'):
            count_cells += 1
            machine = file.split('_')[2]
            try:
                int(machine) # if possible, extract last three digit only
                machine = machine[-3:]
            except:
                pass # use the whole, probably is 'sherri'
            machines.add(machine)
    return {'machines': ",".join(machines), "cells":str(count_cells)}



def extract_ZMW_info(smrtpipe_dir):
    """
    ex: smrtpipe/Christoph_Clontech_1to2k/data/
    look at the # of .ccs.h5 files and assume 150k for now to get the P1
    """
    zmw_per_h5 = 50000 # hardcoded assumption
    files = fnmatch.filter(os.listdir(os.path.join(smrtpipe_dir, 'data')), '*.ccs.h5')
    total_zmw = len(files) * zmw_per_h5
    p1 = int(os.popen("grep -c \">\" " + os.path.join(smrtpipe_dir, 'data', 'reads_of_insert.fasta')).read().strip())
    return {'P1': "{0:.1f}".format(100.*p1/total_zmw), "ZMWs": str(p1)}
    

def read_evaled_summary(filename):
    """
    NOTE: this is now updated for ReadsOfInsert (even tho still called CCS)

    Number of CCS reads: 7318
    Avg. CCS readlength: 2716
    % of artificial chimeras in CCS: 9/6214 (0.1%)
    ----- CCS only ------
    Total number of ZMWs: 6205
    Total number of unmapped: 3 (0.0%)
    Avg. coverage: 99.4%
    Avg. observed accuracy: 95.1%
    Chimera rate: 0.9%
    """
    print >> sys.stderr, "processing", filename
    d = {}
    with open(filename) as f:
        while True:
            if f.tell() == os.stat(filename).st_size: break
            line = f.readline().strip()
            if line.startswith('Number of CCS reads:'):
                d['ccs_yield'] = line.strip().split(': ')[1]
            elif line.startswith('Avg. CCS readlength:'):
                d['ccs_avg_len'] = line.strip().split(': ')[1]
            elif line.startswith('% of artificial chimeras'):
                a, b = line.strip().split(': ')[1].split()
                a, c = a.split('/')
                d['art_chimera'] = "{0} {1}".format(a, b)
                d['fl_zmws'] = c
                d['fl_good'] = str(int(c) - int(a))
            elif line.startswith('----- CCS only ------'):
                assert f.readline().startswith('Total number of ZMWs:')
                d['gmap_ccs_unmapped'] = f.readline().strip().split(': ')[1]
                #assert f.readline().startswith('Total number of unmapped: ')
                d['gmap_ccs_cov'] = f.readline().strip().split(': ')[1]
                d['gmap_ccs_acc'] = f.readline().strip().split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera rate:')
                d['gmap_ccs_chimera'] = line.split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera from missed adapter: ')
                d['gmap_ccs_chimera_from_adapter'] = line.split(': ')[1]
            elif line.startswith('Number of aligned refs:'):
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_refs'] = "{0} {1}".format(a, b)
            elif line.startswith('Number of aligned ZMWs:'):
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_ZMWs'] = "{0} {1}".format(a, b)
            elif line.startswith('Number of well-aligned refs:'):
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_refs_well'] = "{0} {1}".format(a, b)
            elif line.startswith('Number of well-aligned ZMWs:'):
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['blasr_num_ZMWs_well'] = "{0} {1}".format(a, b)
                break # END OF FILE!
    return d

from collections import defaultdict
def read_primer_summary(filename, is_ccs):
    d = defaultdict(lambda: -1)
    with open(filename) as f:
        assert f.readline().startswith("------ 5' primer seen sumary ----")
        if not is_ccs: f.readline()
        d['5seen'] = f.readline().strip().split()[-1][1:-1]
        if not is_ccs: f.readline()
        assert f.readline().startswith("------ 3' primer seen sumary ----")
        if not is_ccs: f.readline()
        d['3seen'] = f.readline().strip().split()[-1][1:-1]
        if not is_ccs: f.readline()
        assert f.readline().startswith("------ 5'&3' primer seen sumary ----")
        f.readline()
        if not is_ccs: f.readline()
        if not is_ccs: f.readline()
        assert f.readline().startswith("------ 5'&3'&polyA primer seen sumary ----")
        if not is_ccs: f.readline()
        d['53Aseen'] = f.readline().strip().split()[-1][1:-1]
        if not is_ccs: f.readline()
    return d


def collect_info_for_name(name):
    d = {'name':name}
    d['ccs_primer'] = read_primer_summary(os.path.join('primer_match', name, 'reads_of_insert.53Aseen_trimmed_changeid.fa.primer_info.txt.summary'), False)
    d.update(read_evaled_summary(os.path.join('primer_match', name, 'evaled_summary.txt')))
    d.update(extract_run_info(os.path.join('runs', name)))
    d.update(extract_ZMW_info(os.path.join('smrtpipe', name)))
    return d

#fields = ['machines', 'cells', 'ZMWs', 'P1', 'ccs_yield', 'ccs_avg_len', 'nonccs_avg_len', '
def write_header(f):
    f.write("Sample,Machine,Cells,SequencingZMWs,Productivity (%),CCS yield,Avg. CCS readlength,,")
    f.write("5' primer seen,3' primer seen,FL %,")
    f.write("# of FL ZMWs,Artificial Chimeras,,")
    f.write("Input # of ZMWs,# unmapped,Avg. GMAP coverage,Avg. GMAP accuracy,GMAP chimera rate,GMAP chimera from RT,,")
    f.write("# of aligned refs,# of aligned ZMWs,# of well-aligned refs,# of well-aligned ZMWs\n")

def write_line(f, d):
    f.write(d['name'] + ',')
    f.write(d['machines'] + ',')
    f.write(d['cells'] + ',')
    f.write(d['ZMWs'] + ',')
    f.write(d['P1'] + ',')
    f.write(d['ccs_yield'] + ',')
    f.write(d['ccs_avg_len'] + ',')

    f.write(',')
    f.write(d['ccs_primer']['5seen'] + ',')
    f.write(d['ccs_primer']['3seen'] + ',')
    f.write(d['ccs_primer']['53Aseen'] + ',')
    f.write(d['fl_zmws'] + ',')
    f.write(d['art_chimera'] + ',')
    f.write(',')

    f.write(d['fl_good'] + ',')
    f.write(d['gmap_ccs_unmapped'] + ',')
    f.write(d['gmap_ccs_cov'] + ',')
    f.write(d['gmap_ccs_acc'] + ',')
    f.write(d['gmap_ccs_chimera'] + ',')
    f.write(d['gmap_ccs_chimera_from_adapter'] + ',')

    f.write(',')
    f.write(d['blasr_num_refs'] + ',')
    f.write(d['blasr_num_ZMWs'] + ',')
    f.write(d['blasr_num_refs_well'] + ',')
    f.write(d['blasr_num_ZMWs_well'] + '\n')

def main():
    import fnmatch
    f = open('cDNA_summary.table.txt','w')
    write_header(f)
    for name in fnmatch.filter(os.listdir('runs'), '[!f]*'):
        if not os.path.isdir(name): continue
        d = collect_info_for_name(name)
        write_line(f, d)
    f.close()

main()
