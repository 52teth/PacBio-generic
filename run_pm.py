#!/usr/bin/env python
import os, sys

"""
Only the runs/ directory have to be set up
Then it does everything needed for primer_match/

(1) run PBBarcode
(2) run barcode_trimmer.py
(3) run count_5seen.py
"""

barcode_trimmer_option = "" # none means 5'-3'
barcode_trimmer_prefix = "53seen_trimmed" # should correspond to the option above

for name in os.listdir('runs/'):
    smrt_d = os.path.join('smrtpipe', name, 'data')
    fa_all1 = os.path.join(smrt_d, 'filtered_subreads.fasta')
    fa_ccs1 = os.path.join(smrt_d, 'filtered_CCS_subreads.fasta')
    if not os.path.exists(fa_ccs1):
        print >> sys.stderr, "{0} does not yet exist. Skipping {1}".format(fa_ccs1, name)
        continue

    pm_d = os.path.join('primer_match', name)
    if not os.path.exists(pm_d): os.makedirs(pm_d)
    fa_primer = os.path.join(pm_d, 'primers.fa')
    if not os.path.exists(fa_primer):
        if os.path.exists('primer_match/primers.fa'):
            cmd = "cp primer_match/primers.fa " + fa_primer
            if os.system(cmd)!=0:
                print >> sys.stderr, "Trouble with {0}.".format(cmd)
                continue
        else:
            print >> sys.stderr, "Need primer file {0}!".format(fa_primer)
            continue

    fa_all2 = os.path.join(pm_d, 'filtered_subreads.fasta')
    fa_ccs2 = os.path.join(pm_d, 'filtered_CCS_subreads.fasta')
    if not os.path.exists(fa_all2) and os.system("ln -s " + os.path.abspath(fa_all1) + " " + fa_all2)!=0: 
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_all1, name)
        continue
    if not os.path.exists(fa_ccs2) and os.system("ln -s " + os.path.abspath(fa_ccs1) + " " + fa_ccs2)!=0:
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_ccs1, name)
        continue

    out_all = os.path.join(pm_d, 'output')
    out_ccs = os.path.join(pm_d, 'outputCCS')
    if not os.path.exists(out_all):
        cmd = "PacBioBarcodeIDCCS.py {0} {1} {2}".format(fa_all2, fa_primer, out_all)
        if os.system(cmd)!=0:
            print >> sys.stderr, "Trouble running {0}.".format(cmd)
    cmd = "barcode_trimmer.py -i {0} -d {1} {2} -o {3}/filtered_subreads.{4}.fa".format(\
                fa_all2, out_all, barcode_trimmer_option, pm_d, barcode_trimmer_prefix)
    if os.system(cmd)!=0:
        print >> sys.stderr, "Trouble running {0}".format(cmd)
    cmd = "count_5seen.py {0}/filtered_subreads.{2}.fa.primer_info.txt {1} > {0}/filtered_subreads.{2}.fa.primer_info.txt.summary".format(pm_d, fa_all2, barcode_trimmer_prefix)
    os.system(cmd)
    if not os.path.exists(out_ccs):
        cmd = "PacBioBarcodeIDCCS.py {0} {1} {2}".format(fa_ccs2, fa_primer, out_ccs)
        if os.system(cmd)!=0:
            print >> sys.stderr, "Trouble running {0}.".format(cmd)
    cmd = "barcode_trimmer.py -i {0} -d {1} {2} -o {3}/filtered_CCS_subreads.{4}.fa".format(\
                fa_ccs2, out_ccs, barcode_trimmer_option, pm_d, barcode_trimmer_prefix)
    if os.system(cmd)!=0:
        print >> sys.stderr, "Trouble running {0}".format(cmd)
    cmd = "count_5seen.py {0}/filtered_CCS_subreads.{2}.fa.primer_info.txt {1} > {0}/filtered_CCS_subreads.{2}.fa.primer_info.txt.summary".format(pm_d, fa_ccs2, barcode_trimmer_prefix)
    os.system(cmd)

    

