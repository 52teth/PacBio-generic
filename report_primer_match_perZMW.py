import os, sys
from csv import reader, DictReader
from collections import defaultdict

def report_primer_match_perZMW(seen53_filename, redux_filename):
    seqid_to_primermatch = {}
    for r in DictReader(open(redux_filename), delimiter='\t'):
        seqid_to_primermatch[r['fastaid'].split()[0]] = r['barcode'] # the .split is for [revcomp]
        
    pm_count = defaultdict(lambda:0) # primer name --> count
    last_key = (None,None)  # (Movie, HoleNumber)
    last_len = 0
    last_pm = None
    for r in DictReader(open(seen53_filename), delimiter='\t'):
        key = r['Movie'], r['HoleNumber']
        cur_len = int(r['rEnd']) - int(r['rStart'])
        if key != last_key:
            # output from last
            if last_pm is not None: pm_count[last_pm] += 1
            last_len = cur_len
            last_key = key
        elif cur_len > last_len:
            last_len = cur_len
            seqid = r['Movie'] + '/' + r['HoleNumber'] + '/' + str(r['rStart']) + '_' + str(r['rEnd'])
            last_pm = seqid_to_primermatch[seqid]
       
    keys = pm_count.keys()
    keys.sort()
    for k in keys:
        print("{0}:{1}".format(k, pm_count[k]))     
        
if __name__ == "__main__":
    report_primer_match_perZMW(sys.argv[1], sys.argv[2])