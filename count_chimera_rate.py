import os, sys
from Bio import SeqIO

"""
Count chimera rate per ZMW
"""
def count_chimera_rate(is_chimera_filename, non_chimera_filename):
    is_chimera_zmw = set()
    non_chimera_zmw = set()
    for r in SeqIO.parse(open(is_chimera_filename), 'fasta'):
        is_chimera_zmw.add(r.id[:r.id.rfind('/')])
    for r in SeqIO.parse(open(non_chimera_filename), 'fasta'):
        non_chimera_zmw.add(r.id[:r.id.rfind('/')])

    a = len(is_chimera_zmw)
    b = len(non_chimera_zmw.difference(is_chimera_zmw))
    print "Confusing ZMWs (some chimeric some not): {0}".format(len(is_chimera_zmw.intersection(non_chimera_zmw)))
    print "Chimera rate (per ZMW): {0}/{1} {2:.2f}".format(a, a+b, 100.*a/(a+b))
    return is_chimera_zmw, non_chimera_zmw


if __name__ == "__main__":
    for d in os.listdir('.'):
        if os.path.isdir(d):
            print "Directory", d
            f1 = os.path.join(d,'filtered_subreads.53seen_trimmed.fa.is_chimera.fa')
            f2 = os.path.join(d,'filtered_subreads.53seen_trimmed.fa.non_chimera.fa')
            if os.path.exists(f1):
                count_chimera_rate(os.path.join(d,'filtered_subreads.53seen_trimmed.fa.is_chimera.fa'),\
                        os.path.join(d,'filtered_subreads.53seen_trimmed.fa.non_chimera.fa'))
