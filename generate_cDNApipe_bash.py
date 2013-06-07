import os, sys


def get_nonCCS(bash_f):
    bash_f.write("grab_nonCCS_subreads.py filtered_subreads.fasta ccs_reads.fasta nonccs_subreads.fasta > ccs_to_nonccs.txt || exit $?\n")
    
def get_primer_info(bash_f, args):
    """
    params: --change-seqid, --must-see-polyA
    """
    bash_f.write("hmmer_wrapper.py -d output -i nonccs_subreads.fasta --change-seqid --must-see-polyA -o nonccs_subreads.53Aseen_trimmed_changeid.fa --cpus {cpus} || exit $?\n".format(cpus=args.cpus))
    bash_f.write("summarize_primer_info.py nonccs_subreads.53Aseen_trimmed_changeid.fa.primer_info.txt > nonccs_subreads.53Aseen_trimmed_changeid.fa.primer_info.txt.summary\n")
    
    bash_f.write("hmmer_wrapper.py -d outputCCS -i ccs_reads.fasta --change-seqid --must-see-polyA -o ccs_reads.53Aseen_trimmed_changeid.fa --cpus {cpus} || exit $?\n".format(cpus=args.cpus))
    bash_f.write("summarize_primer_info.py ccs_reads.53Aseen_trimmed_changeid.fa.primer_info.txt > ccs_reads.53Aseen_trimmed_changeid.fa.primer_info.txt.summary\n")
    
def get_chimera(bash_f, args):
    bash_f.write("chimera_finder.py -d outputNonCCS_chimera --cpus {cpus} -i nonccs_subreads.53Aseen_trimmed_changeid.fa || exit $?\n".format(cpus=args.cpus))
    bash_f.write("chimera_finder.py -d outputCCS_chimera --cpus {cpus} -i ccs_reads.53Aseen_trimmed_changeid.fa || exit $?\n".format(cpus=args.cpus))
        
def run_BLASR(bash_f, args):
    bash_f.write("pbalign.py --maxHits 10 --hitPolicy all --algorithmOptions \"-nproc {cpus} -bestn 10 -nCandidates 10\" nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa {ref} nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa.blasr.sam\n".format(cpus=args.cpus, ref=args.ref))
    bash_f.write("pbalign.py --maxHits 10 --hitPolicy all --algorithmOptions \"-nproc {cpus} -bestn 10 -nCandidates 10\"  ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa  {ref} ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa.blasr.sam\n".format(cpus=args.cpus, ref=args.ref))    

def run_gmap(bash_f, args):
    bash_f.write("gmap -D ~/share/gmap_db/ -d {gmap_db} -t {cpus} -f gff3_gene -n 0 nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa > nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff 2> nonccs_subreads.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff.log\n".format(cpus=args.cpus, gmap_db=args.gmap_db))
    bash_f.write("gmap -D ~/share/gmap_db/ -d {gmap_db} -t {cpus} -f gff3_gene -n 0 ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa > ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff 2> ccs_reads.53Aseen_trimmed_changeid.fa.non_chimera.fa.gff.log\n".format(cpus=args.cpus, gmap_db=args.gmap_db))
    

def main(args):
    f = open('cDNApipe.sh', 'w')
    f.write("#!/bin/bash\n")
    f.write(". /home/UNIXHOME/etseng/.VENV/bin/activate\n")
    
    get_nonCCS(f)
    get_primer_info(f, args)
    get_chimera(f, args)
    run_BLASR(f, args)
    run_gmap(f, args)
    
    f.close()       
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--ref", default='/mnt/secondary/Share/Smrtanalysis-alpha/opt/smrtanalysis/common/references/Gencode15/', help="Transcript reference location (default: Gencode15)")
    parser.add_argument("--gmap_db", default='hg19', help="GMAP db name (default: hg19)")
    parser.add_argument("--cpus", type=int, default=24, help="Number of CPUs (default: 24)")
    
    args = parser.parse_args()
    main(args)     