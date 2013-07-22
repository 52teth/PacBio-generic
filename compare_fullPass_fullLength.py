import os, sys
from csv import DictReader, DictWriter
from pbcore.io import BasH5Reader

def get_fullpass_info(zmw):
    """
    Given a zmw object from BasH5Reader
    Return the list of subreads that are full-pass
    """
    if len(zmw.adapters) == 0: # there are no adapters
        return [], (-1,-1), (-1, -1)

    adap = []
    for a, b in zmw.adapterRegions:
        adap.append(a)
        adap.append(b)
    n = len(adap)

    i = 0
    result = []
    first_pass = (-1, -1)
    last_pass = (-1, -1)
    if zmw.subreads[0].readStart < adap[0]: 
        first_pass = (zmw.subreads[0].readStart, zmw.subreads[0].readEnd)
    if len(zmw.subreads) > 1 and zmw.subreads[-1].readStart >= adap[-1]:
        last_pass = (zmw.subreads[-1].readStart, zmw.subreads[-1].readEnd)
    for x in zmw.subreads:
        try:
            i = adap.index(x.readStart, i)
            if i >= 0 and i < n - 1 and adap[i+1] == x.readEnd:
                result.append((x.readStart, x.readEnd))
        except:
            pass
    return result, first_pass, last_pass

def uneven_subread_lengths(fp_list):
    """
    Return True if the subread lengths from this ZMW are uneven
    (indicating missed or mis-called adapters)
    """
    import numpy as np
    if len(fp_list) <= 1: return False
    lens = np.array([x[1]-x[0] for x in fp_list])
    return max(lens) > 1.5* min(lens)

def is_full_pass(fp_list, s, e):
    """
    fp_list -- must be sorted
    """
    for start, end in fp_list:
        if start <= s < e <= end: return True
    return False

def is_first_pass(fp_list, s, e):
    if len(fp_list) == 0: return False
    start, end = fp_list[0]
    if start <= s < e <= end: return True

def enhance_primer_info(primer_info_filename, bas_dict):
    fields = ['ID', 'strand', '5seen', 'polyAseen', '3seen', '5end', 'polyAend', '3end', 'primer', 'firstPass', 'lastPass', 'fullPass', 'unevenZMW']
    fields1 = ['ZMW', 'has_FL', 'has_FP', 'FL_FP_discordant_exclude1PLP']
    f = open(primer_info_filename + '.enhanced.txt', 'w')
    f1 = open(primer_info_filename + '.enhanced_byZMW.txt', 'w')
    f.write("\t".join(fields) + '\n')
    f1.write("\t".join(fields1) + '\n')
    writer = DictWriter(f, fields, delimiter='\t')
    last_movie = None
    last_hole = None
    last_fp_list = None
    last_uneven = None
    last_has_FP = None
    last_has_FL = None
    last_FL_FP_bad = None
    for r in DictReader(open(primer_info_filename), delimiter='\t'):
        movie, hole, s_e = r['ID'].split('/')
        s, e = map(int, s_e.split('_'))
        if s > e: s, e = e, s
        if movie != last_movie or last_hole != hole:
            if last_has_FP is not None:
                f1.write("{0}\t{1}\t{2}\t{3}\n".format(last_movie+'/'+last_hole, '1' if last_has_FL else '0', '1' if last_has_FP else '0', '1' if last_FL_FP_bad else '0'))
            last_movie = movie
            last_hole = hole
            last_fp_list, last_first_pass, last_last_pass = get_fullpass_info(bas_dict[movie][int(hole)])
            last_uneven = uneven_subread_lengths(last_fp_list)
            last_has_FL = False
            last_has_FP = False
            last_FL_FP_bad = False
        is_fp = is_full_pass(last_fp_list, s, e)
        r['firstPass'] = '1' if last_first_pass[0] <= s < e <= last_first_pass[1] else '0'
        r['lastPass'] = '1' if last_last_pass[0] <= s < e <= last_last_pass[1] else '0'
        r['fullPass'] = '1' if is_fp else '0'
        r['unevenZMW'] = '1' if last_uneven else '0'
        is_fl = r['5seen']=='1' and r['3seen']=='1'
        if is_fl: last_has_FL = True
        if is_fp: last_has_FP = True
        if (is_fl ^ is_fp) and r['firstPass']=='0' and r['lastPass']=='0': last_FL_FP_bad = True
        writer.writerow(r)
    f.close()
    f1.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("primer_info_filename")
    parser.add_argument("input_fofn")
    args = parser.parse_args()

    bas_dict = {}
    for line in open(args.input_fofn):
        file = line.strip()
        if file.endswith('.bas.h5'): 
            movie = os.path.basename(file)[:-7]
        else: 
            movie = os.path.basename(file)[:-9]  #.1.bax.h5
            file  = file[:-9] + '.bas.h5'
        if movie not in bas_dict:
            print >> sys.stderr, "Reading {0}".format(file)
            bas_dict[movie] = BasH5Reader(file)

    enhance_primer_info(args.primer_info_filename, bas_dict)

