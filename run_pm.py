#!/usr/bin/env python
import os, sys, subprocess

"""
Only the runs/ directory have to be set up
Then it does everything needed for primer_match/

(1) run hmmer_wrapper.py
(2) run summarize_primer_info.py
"""
species = sys.argv[1]

barcode_trimmer_option = "--change-seqid, --must-see-polyA --min-seqlen 100 --cpus 12" # none means 5'-3'
barcode_trimmer_prefix = "53Aseen_trimmed" # should correspond to the option above

if species == 'hg19':
    cmd_cDNA = "generate_cDNApipe_bash.py --cpus 12 --cmd_filename cDNApipe.sh"
else:
    cmd_cDNA = "generate_cDNApipe_bash.py --ref /mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/rat_UCSC --gmap_db rn5 --cpus 12 --cmd_filename cDNApipe.sh"

for name in os.listdir('runs/'):
    smrt_d = os.path.join('smrtpipe', name, 'data')
    fa_ccs1 = os.path.join(smrt_d, 'reads_of_insert.fasta')
    if not os.path.exists(fa_ccs1):
        print >> sys.stderr, "{0} does not yet exist. Skipping {1}".format(fa_ccs1, name)
        continue

    pm_d = os.path.join('primer_match', name)
    if not os.path.exists(pm_d): os.makedirs(pm_d)
    else: continue # skipping pm_d because already exists
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

    fa_ccs2 = os.path.join(pm_d, 'reads_of_insert.fasta')
    if not os.path.exists(fa_ccs2) and os.system("ln -s " + os.path.abspath(fa_ccs1) + " " + fa_ccs2)!=0:
        print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(fa_ccs1, name)
        continue


    cwd = os.popen("pwd").read().strip()
    os.chdir(pm_d)
    subprocess.check_call(cmd_cDNA, shell=True)
    cmd = "qsub -cwd -S /bin/bash -pe smp 12 cDNApipe.sh"
    print >> sys.stderr, "submitting job for ", pm_d
    subprocess.check_call(cmd, shell=True)
    os.chdir(cwd)
