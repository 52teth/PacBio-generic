#!/usr/bin/env python
"""
Sets up directories in smrtpipe/ according to runs/ directory

CONFIG must have:
 --reads_of_insert_protocol
 --rs_filter_protocol
 --cpus
"""
import os, sys, subprocess

CONFIG_PARAMS = ['reads_of_insert_protocol', 'rs_filter_protocol', 'cpus', 'gmap_db', 'ref_transcript', 'primer_filename']

def run_cmd(cmd):
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "problem running", cmd, ". Abort!"
        sys.exit(-1)
        
def read_config(filename):
    config = {}
    for line in open(filename):
        if not line.startswith('#') and line.find('=') > 0:
            k, v = line.strip().split('=')
            config[k.strip()] = v.strip()
            
    for k in CONFIG_PARAMS:
        if k not in config:
            print >> sys.stderr, "Parameter {0} not found in {1}. Abort!".format(k, filename)
            sys.exit(-1)
            
    return config
            
        
def setup_smrtpipe(name, config, protocol='reads_of_insert'):
    """
    """
    d_run = os.path.join('runs', name)
    if not os.path.isdir(d_run):
        print >> sys.stderr, "Directory {0} does not exist. Abort!".format(d_run)
        sys.exit(-1)
    
    d = os.path.join('smrtpipe/', name, protocol)
    if not os.path.exists(d):
        print >> sys.stderr, "creating directory", d
        os.makedirs(d)
    cmd = "find {0}/*/Analysis_Results/*.bax.h5 > {1}/input.fofn".format(os.path.abspath(d_run),d)
    run_cmd(cmd)
        
    cmd = "fofnToSmrtpipeInput.py {0}/input.fofn > {0}/input.xml".format(d)
    run_cmd(cmd)
    
    cmd = "cp {0} {1}/settings.xml".format(config[protocol+'_protocol'], d)
    run_cmd(cmd)

    cmd = "ln -s {0}/common {1}/common".format(os.environ['SEYMOUR_HOME'], d)
    run_cmd(cmd)
    
    with open(os.path.join(d, 'smrtpipe.sh'), 'w') as f:
        f.write("#/!bin/bash\n")
        f.write("smrtpipe.py --params={0}/settings.xml xml:{0}/input.xml -D NPROC={1} --output {0}\n".format(os.path.abspath(d), config['cpus']))
        
    return os.path.abspath(f.name)

def setup_primer_match(name, config):
    d = os.path.join('primer_match', name)
    if not os.path.exists(d):
        print >> sys.stderr, "creating directory", d
        os.makedirs(d)
        
    cmd = "cp {0} {1}".format(config['primer_filename'], d)
    run_cmd(cmd)  
    return d  
        
def main(name):
    config = read_config('cDNA_smrtpipe.config')
    fname1 = setup_smrtpipe(name, config, 'reads_of_insert')
    fname2 = setup_smrtpipe(name, config, 'rs_filter')
    
    pmdir = setup_primer_match(name, config)
    
    print "bash", fname1
    print "bash", fname2
    print "cd", pmdir
    print "ln -s {0}/data/reads_of_insert.fasta reads_of_insert.fasta".format(os.path.dirname(fname1))
    print "ln -s {0}/data/filtered_subreads.fasta filtered_subreads.fasta".format(os.path.dirname(fname2))
        
    pm_cmd = 'cDNApipe.sh'
    print "generate_cDNApipe_bash.py --ref {r} --gmap_db {g} --cpus {c} --cmd_filename {o}".format(\
        r=config['ref_transcript'], g=config['gmap_db'], c=config['cpus'], o=pm_cmd) 
    
    print "bash", pm_cmd

    
if __name__ == "__main__": 
    name = sys.argv[1]
    main(name)    
    
