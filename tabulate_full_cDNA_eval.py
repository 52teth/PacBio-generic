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
    ex: smrtpipe/Christoph_Clontech_1to2k/data/filtered_summary.csv
    Movie,ReadId,#Bases,Readlength,ReadScore,SequencingZMW,Productivity,PassedFilter
    """
    all = 0
    p1 = 0
    file = os.path.join(smrtpipe_dir, 'data/filtered_summary.csv')
    for r in DictReader(open(file), delimiter=','):
        all += 1
        if r['SequencingZMW']=='1' and r['Productivity']=='1' and r['PassedFilter']=='1':
            p1 += 1

    return {'P1': "{0:.1f}".format(100.*p1/all), "ZMWs": str(p1)}

def read_evaled_summary(filename):
    """
    % of ZMWs with CCS: 7318/54122 (13.5%)
    Number of CCS reads: 7318
    Number of non-CCS subreads: 71230
    Avg. CCS readlength: 2716
    Avg. non-CCS subread readlength: 1637
    % of artificial chimeras in CCS: 9/6214 (0.1%)
    % of artificial chimeras in non-CCS subreads: 2/2904 (0.1%)
    % of artificial chimeras: 11/9118 (0.1%)
    ----- CCS only ------
    Total number of ZMWs: 6205
    Total number of unmapped: 3 (0.0%)
    Avg. coverage: 99.4%
    Avg. observed accuracy: 95.1%
    Chimera rate: 0.9%
    ----- non-CCS subread only ------
    Total number of ZMWs: 2902
    Total number of unmapped: 157 (5.4%)
    Avg. coverage: 69.3%
    Avg. observed accuracy: 82.0%
    Chimera rate: 3.1%
    ----- Combined -------
    Total number of ZMWs: 9107
    Total number of unmapped: 160 (1.8%)
    Avg. coverage: 89.8%
    Avg. observed accuracy: 90.9%
    Chimera rate: 1.6%
    """
    d = {}
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            if line.startswith('% of ZMWs with CCS:'):
                a, b = line.strip().split(': ')[1].split()
                a = a.split('/')[0]
                d['ccs_yield'] = "{0} {1}".format(a, b)
            elif line.startswith('Avg. CCS readlength:'):
                d['ccs_avg_len'] = line.strip().split(': ')[1]
            elif line.startswith('Avg. non-CCS subread readlength:'):
                d['nonccs_avg_len'] = line.strip().split(': ')[1]
            elif line.startswith('% of artificial chimeras:'):
                a, b = line.strip().split(': ')[1].split()
                a, c = a.split('/')
                d['art_chimera'] = "{0} {1}".format(a, b)
                d['fl_zmws'] = c
            elif line.startswith('----- CCS only ------'):
                assert f.readline().startswith('Total number of ZMWs:')
                assert f.readline().startswith('Total number of unmapped: ')
                assert f.readline().startswith('Avg. coverage:')
                assert f.readline().startswith('Avg. observed accuracy:')
                line = f.readline().strip()
                assert line.startswith('Chimera rate:')
                d['gmap_ccs_chimera'] = line.split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera from missed adapter: ')
                d['gmap_ccs_chimera_from_adapter'] = line.split(': ')[1]
            elif line.startswith('----- non-CCS subread only ------'):
                assert f.readline().startswith('Total number of ZMWs:')
                assert f.readline().startswith('Total number of unmapped: ')
                assert f.readline().startswith('Avg. coverage:')
                assert f.readline().startswith('Avg. observed accuracy:')
                line = f.readline().strip()
                assert line.startswith('Chimera rate:')
                d['gmap_nonccs_chimera'] = line.split(': ')[1]
                line = f.readline().strip()
                assert line.startswith('Chimera from missed adapter: ')
                d['gmap_nonccs_chimera_from_adapter'] = line.split(': ')[1]
            elif line.startswith('----- Combined -------'):
                d['fl_nonchimera_zmws'] = f.readline().strip().split(': ')[1] # Total number of ZMWs:
                d['gmap_unmapped'] = f.readline().strip().split(': ')[1] # Total number of unmapped: 160 (1.8%)
                d['gmap_avg_cov'] = f.readline().strip().split(': ')[1] # Avg. coverage: 89.8%
                d['gmap_avg_acc'] = f.readline().strip().split(': ')[1] # Avg. observed accuracy: 90.9%
                d['gmap_chimera'] = f.readline().strip().split(': ')[1] # Chimera rate: 1.6%
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


def read_primer_summary(filename, is_ccs):
    d = {}
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
    d['ccs_primer'] = read_primer_summary(os.path.join('primer_match', name, 'ccs_reads.53Aseen_trimmed_changeid.fa.primer_info.txt.summary'), True)
    d['nonccs_primer'] = read_primer_summary(os.path.join('primer_match', name, 'nonccs_subreads.53Aseen_trimmed_changeid.fa.primer_info.txt.summary'), False)
    d.update(read_evaled_summary(os.path.join('primer_match', name, 'evaled_summary.txt')))
    d.update(extract_run_info(os.path.join('runs', name)))
    d.update(extract_ZMW_info(os.path.join('smrtpipe', name)))
    return d

#fields = ['machines', 'cells', 'ZMWs', 'P1', 'ccs_yield', 'ccs_avg_len', 'nonccs_avg_len', '
def write_header(f):
    f.write("Sample,Machine,Cells,SequencingZMWs,Productivity (%),CCS yield,Avg. CCS readlength,Avg. nonCCS readlength,,CCS 5' primer seen,CCS 3' primer seen,CCS FL %,non-CCS 5' primer seen,non-CCS 3' primer seen,non-CCS FL%,# of FL ZMWs,Artificial Chimeras,,Input # of ZMWs,# unmapped,Avg. GMAP coverage,Avg. GMAP accuracy,GMAP chimera rate,GMAP CCS chimera,GMAP CCS chimera from adapter,GMAP non-CCS chimera,GMAP non-CCS chimera from adapter,,# of aligned refs,# of aligned ZMWs,# of well-aligned refs,# of well-aligned ZMWs\n")

def write_line(f, d):
    f.write(d['name'] + ',')
    f.write(d['machines'] + ',')
    f.write(d['cells'] + ',')
    f.write(d['ZMWs'] + ',')
    f.write(d['P1'] + ',')
    f.write(d['ccs_yield'] + ',')
    f.write(d['ccs_avg_len'] + ',')
    f.write(d['nonccs_avg_len'] + ',')
    f.write(',')
    f.write(d['ccs_primer']['5seen'] + ',')
    f.write(d['ccs_primer']['3seen'] + ',')
    f.write(d['ccs_primer']['53Aseen'] + ',')
    f.write(d['nonccs_primer']['5seen'] + ',')
    f.write(d['nonccs_primer']['3seen'] + ',')
    f.write(d['nonccs_primer']['53Aseen'] + ',')
    f.write(d['fl_zmws'] + ',')
    f.write(d['art_chimera'] + ',')
    f.write(',')
    f.write(d['fl_nonchimera_zmws'] + ',')
    f.write(d['gmap_unmapped'] + ',')
    f.write(d['gmap_avg_cov'] + ',')
    f.write(d['gmap_avg_acc'] + ',')
    f.write(d['gmap_chimera'] + ',')
    f.write(d['gmap_ccs_chimera'] + ',')
    f.write(d['gmap_ccs_chimera_from_adapter'] + ',')
    f.write(d['gmap_nonccs_chimera'] + ',')
    f.write(d['gmap_nonccs_chimera_from_adapter'] + ',')
    f.write(',')
    f.write(d['blasr_num_refs'] + ',')
    f.write(d['blasr_num_ZMWs'] + ',')
    f.write(d['blasr_num_refs_well'] + ',')
    f.write(d['blasr_num_ZMWs_well'] + '\n')

def main():
    f = open('cDNA_summary.table.txt','w')
    write_header(f)
    for name in os.listdir('runs'):
        d = collect_info_for_name(name)
        write_line(f, d)
    f.close()

main()
