#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Collect .ccs.bam into ccs.fasta/fastq
Run pbtranscript.py classify
"""

import os, sys, glob

def check_file_or_run(filename):
    def func(f):
        def g(bash_f, args):
            if args.force or not os.path.exists(filename):
                f(bash_f, args)
            else:
                print >> sys.stderr, "{0} already exists. No need to run {1}.".format(filename, f.__name__)
        return g
    return func

@check_file_or_run('ccs.fasta')
def collect_ccs_bam(bash_f, args):
    for file in glob.glob('*.ccs.bam'):
        movie = file.split('.')[0]
        assert os.path.exists(movie + '.ccs.report')
        bash_f.write("bamtools convert -in {0}.ccs.bam -format fastq > {0}.ccs.fastq\n".format(movie))
    bash_f.write("cat *.ccs.fastq > ccs.fastq\n")
    bash_f.write("fq2fa.py ccs.fastq\n")

@check_file_or_run('isoseq_flnc.fasta')
def generate_classify(bash_f, args):
    """
    Run classify
    Convert to FASTQ files
    """
    cmd = "pbtranscript.py classify ccs.fasta isoseq_draft.fasta --flnc isoseq_flnc.fasta --nfl isoseq_nfl.fasta " + \
        "--cpus {c} -p {p} --min_seq_len {m}\n".format(c=args.cpus, p=args.primer, m=args.min_seq_len)
    bash_f.write(cmd)
    bash_f.write("python ~/GitHub/PB_ICE/get_fq_from_fa.py_log\n")


def main(args):
    f = open(args.cmd_filename, 'w')
    f.write("#!/bin/bash\n")
    f.write(". {0}/bin/activate\n".format(os.environ['VIRTUAL_ENV']))

    collect_ccs_bam(f, args)
    generate_classify(f, args)

    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--min_seq_len", type=int, default=300, help="Minimum seq length (default: 300 bp)")
    parser.add_argument("-p", "--primer", default='/home/UNIXHOME/etseng/GitHub/cDNA_primer/pbtranscript-tofu//pbtranscript//pbtools/pbtranscript//data/primers.fa', help="Primer filename")
    parser.add_argument("--cpus", type=int, default=12, help="Number of CPUs (default: 12)")
    parser.add_argument("--force", action="store_true", default=False, help="Run all commands even if output files are present for some steps")
    parser.add_argument("--cmd_filename", default='ccs_classify.sh', help="Output CMD filename")

    args = parser.parse_args()
    main(args)
