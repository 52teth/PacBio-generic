#!/usr/bin/env python
import os, sys

"""
Only the runs/ directory have to be set up
Then it does everything needed for primer_match/

(1) run hmmer_wrapper.py
(2) run summarize_primer_info.py
"""

barcode_trimmer_option = "--cpus 24" # none means 5'-3'
barcode_trimmer_prefix = "53seen_trimmed" # should correspond to the option above

for name in os.listdir('runs/'):
    smrt_d = os.path.join('smrtpipe', name, 'data')
    fa_all1 = os.path.join(smrt_d, 'filtered_subreads.fasta')
    fa_ccs1 = os.path.join(smrt_d, 'ccs_reads.fasta')
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
    fa_ccs2 = os.path.join(pm_d, 'ccs_reads.fasta')
    if not os.path.exists(fa_all2) and os.system("ln -s " + os.path.abspath(fa_all1) + " " + fa_all2)!=0: 
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_all1, name)
        continue
    if not os.path.exists(fa_ccs2) and os.system("ln -s " + os.path.abspath(fa_ccs1) + " " + fa_ccs2)!=0:
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_ccs1, name)
        continue

    out_all = os.path.join(pm_d, 'output')
    out_ccs = os.path.join(pm_d, 'outputCCS')
    pf = os.path.join(pm_d, 'primers.fa')
    if not os.path.exists(out_all):
        output_filename = os.path.join(pm_d, "filtered_subreads.{op}.fa".format(op=barcode_trimmer_prefix))
        cmd = "hmmer_wrapper.py -d {d} -i {i} -p {pf} -o {o} {e}".format(d=out_all, i=fa_all2, o=output_filename, pf=pf, e=barcode_trimmer_option)
        if os.system(cmd)!=0:
            print >> sys.stderr, "Trouble running {0}.".format(cmd)
        cmd = "summarize_primer_info.py {0}.primer_info.txt > {0}.primer_info.txt.summary".format(output_filename)
        os.system(cmd)
    if not os.path.exists(out_ccs):
        output_filename = os.path.join(pm_d, "ccs_reads.{op}.fa".format(op=barcode_trimmer_prefix))
        cmd = "hmmer_wrapper.py -d {d} -i {i} -p {pf} -o {o} {e}".format(d=out_ccs, i=fa_ccs2, o=output_filename, pf=pf, e=barcode_trimmer_option)
        if os.system(cmd)!=0:
            print >> sys.stderr, "Trouble running {0}.".format(cmd)
        cmd = "summarize_primer_info.py {0}.primer_info.txt > {0}.primer_info.txt.summary".format(output_filename)
        os.system(cmd)

    

