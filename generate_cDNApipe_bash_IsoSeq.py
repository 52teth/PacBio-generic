#!/usr/bin/env python
import os, sys

def check_file_or_run(filename):
    def func(f):
        def g(bash_f, args):
            if args.force or not os.path.exists(filename):
                f(bash_f, args)
            else:
                print >> sys.stderr, "{0} already exists. No need to run {1}.".format(filename, f.__name__)
        return g
    return func

@check_file_or_run('isoseq_flnc.fasta.blasr.sam')
def run_BLASR(bash_f, args):
    if args.ref != 'NA':
        bash_f.write("pbalign.py --maxHits 10 --hitPolicy all --algorithmOptions \"-nproc {cpus} -bestn 10 -nCandidates 10\" isoseq_flnc.fasta {ref} isoseq_flnc.fasta.blasr.sam\n".format(cpus=args.cpus, ref=args.ref))    

@check_file_or_run('isoseq_flnc.fasta.gff.log')
def run_gmap(bash_f, args):
    if args.gmap_db != 'NA':
        bash_f.write("~/bin/gmap -D {gmap_db_dir} -d {gmap_db} -t {cpus} -f samse -n 0 isoseq_flnc.fasta > isoseq_flnc.fasta.sam 2> isoseq_flnc.fasta.sam.log\n".format(cpus=args.cpus, gmap_db=args.gmap_db, gmap_db_dir=args.gmap_db_dir))
        bash_f.write("sort -k 3,3 -k 4,4n isoseq_flnc.fasta.sam > isoseq_flnc.fasta.sorted.sam\n")
        bash_f.write("collapse_isoforms_by_sam.py --input isoseq_flnc.fasta -s isoseq_flnc.fasta.sorted.sam -c 0.99 -i 0.95 -o isoseq_flnc.5merge_c99i95\n")
    

def main(args):
    f = open(args.cmd_filename, 'w')
    f.write("#!/bin/bash\n")
    f.write(". {0}/bin/activate\n".format(os.environ['VIRTUAL_ENV']))
    
    run_BLASR(f, args)
    run_gmap(f, args)
    
    f.close()       
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--ref", default='/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/', help="Transcript reference location (default: Gencode15)")
    parser.add_argument("--gmap_db", default='hg19', help="GMAP db name (default: hg19)")
    parser.add_argument("--gmap_db_dir", default='/home/UNIXHOME/etseng/share/gmap_db/', help="GMAP db directory")
    parser.add_argument("--cpus", type=int, default=24, help="Number of CPUs (default: 24)")
    parser.add_argument("--force", action="store_true", default=False, help="Run all commands even if output files are present for some steps")
    parser.add_argument("--cmd_filename", default='cDNApipe.sh', help="Output CMD filename")
    
    args = parser.parse_args()
    main(args)     
