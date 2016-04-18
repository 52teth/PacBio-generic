#!/usr/bin/env python
import os, sys
import subprocess

def run_cmd(cmd):
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "Error cmd:", cmd



input = sys.argv[1] # ex: test.sam

assert input.endswith('.sam')
prefix = input[:-4]

run_cmd("samtools view -bS {0}.sam > {0}.bam".format(prefix))
run_cmd("samtools sort {0}.bam {0}".format(prefix))
run_cmd("samtools index {0}.bam".format(prefix))


