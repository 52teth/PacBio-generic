#!/usr/bin/env python
import os, sys, subprocess

species = sys.argv[1]

# ----------------------- SETTINGS --------------------- #
NUM_CPUS = 12
gmap_db_dir = '/home/UNIXHOME/etseng/share/gmap_db_new'

# option for GMAP and BLASR
if species == 'hg19':
    transcript_ref_dir = '/home/UNIXHOME/etseng/share/gencode/'
    transcript_ref_fasta = '/home/UNIXHOME/etseng/share/gencode/gencode.v23.pc_transcripts.derep.fa'
    gmap_db_name = 'hg19'
elif species == 'hg19_gmaponly':
    transcript_ref_dir = 'NA'
    transcript_ref_fasta = 'NA'
    gmap_db_name = 'hg19'
elif species == 'chimp':
    transcript_ref_dir = '/home/UNIXHOME/etseng/share/gencode/'
    transcript_ref_fasta = '/home/UNIXHOME/etseng/share/gencode/gencode.v23.pc_transcripts.derep.fa'    
    gmap_db_name = 'panTro4'
elif species == 'rn5':
    transcript_ref_dir = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/rat_UCSC'
    transcript_ref_fasta = '/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/rat_UCSC/sequence/rat_UCSC.fasta'
    gmap_db_name = 'rn5'
elif species == 'mm10':
    transcript_ref_dir = '/pbi/dept/bifx/etseng/genomes/RefSeq'
    transcript_ref_fasta = '/pbi/dept/bifx/etseng/genomes/RefSeq/mouse_RefSeq.rna.fasta'
    gmap_db_name = 'mm10'
elif species == 'neurospora':
    transcript_ref_dir = 'NA'
    transcript_ref_fasta = 'NA'
    gmap_db_name = 'neurospora'
elif species == 'skip':
    transcript_ref_dir = 'NA'
    transcript_ref_fasta = 'NA'
    gmap_db_name = 'NA'
else:
    print >> sys.stderr, "species not specified or unknown! quit!"
    sys.exit(-1)

cmd_cDNA = "generate_cDNApipe_bash_IsoSeq.py --ref {0} --gmap_db {1} --gmap_db_dir {2} --cpus {3} --cmd_filename cDNApipe.sh".format(transcript_ref_fasta, gmap_db_name, gmap_db_dir, NUM_CPUS)
if species == 'skip':
    cmd_eval = "eval_cDNApipe_results_IsoSeq.py --skip-GMAP --skip-BLASR > evaled_summary.txt"
else:
    cmd_eval = "eval_cDNApipe_results_IsoSeq.py --ref_fasta_filename {0} > evaled_summary.txt".format(transcript_ref_fasta)

# ----------------------- SETTINGS --------------------- #

# SANITY CHECK for runs/, smrtpipe/, eval
if not os.path.exists('runs/'):
    print >> sys.stderr, "runs/ directory does not exist. Quit!"
    sys.exit(-1)
if not os.path.exists('smrtpipe/'):
    print >> sys.stderr, "smrtpipe/ directory does not exist. Quit!"
    sys.exit(-1)
if not os.path.exists('eval/'):
    print >> sys.stderr, "eval/ directory does not exist. Quit!"
    sys.exit(-1)

# FOR EACH <name> in runs/
#  find the corressponding smrtpipe/<name>
#  symbolically link isoseq_flnc.fasta to eval/<name>
#  generate cDNA script
#  submit cDNA script
for name in os.listdir('runs/'):
    smrt_d = os.path.join('smrtpipe', name, 'data')
    files = ['isoseq_flnc.fasta', 'isoseq_nfl.fasta', 'classify_summary.txt', 'isoseq_primer_info.csv', 'reads_of_insert.fasta']
    if not os.path.exists(os.path.join(smrt_d, files[0])):
        print >> sys.stderr, "{0} does not yet exist. Skipping {1}".format(files[0], name)
        continue
    
    pm_d = os.path.join('eval', name)
    if not os.path.exists(pm_d): os.makedirs(pm_d)
    else: continue # skipping pm_d because already exists

    for file in files:
        file_src = os.path.join(smrt_d, file)
        file_path = os.path.join(pm_d, file)
        if not os.path.exists(file_path) and os.system("ln -s " + os.path.abspath(file_src) + " " + file_path)!=0:
            print >> sys.stderr, "Trouble linking {0}. Skipping {1}".format(file_path, name)
            continue
    
    smrt_d2 = os.path.join('smrtpipe', name, 'results')
    files2 = ['fulllength_nonchimeric_readlength_hist.png','roi_npasses_hist.png',\
            'roi_accuracy_hist.png','roi_readlength_hist.png']
    for file in files2:
        file_src = os.path.join(smrt_d2, file)
        file_path = os.path.join(pm_d, name+'.'+file)
        os.system("ln -s " + os.path.abspath(file_src) + " " + file_path)


    cwd = os.popen("pwd").read().strip()
    os.chdir(pm_d)
    subprocess.check_call(cmd_cDNA, shell=True)
    with open('cDNApipe.sh', 'a') as f: 
        f.write(cmd_eval + '\n')
    cmd = "qsub -cwd -V -S /bin/bash -pe smp {cpus} cDNApipe.sh".format(cpus=NUM_CPUS)
    print >> sys.stderr, "submitting job for ", pm_d
    subprocess.check_call(cmd, shell=True)
    os.chdir(cwd)
