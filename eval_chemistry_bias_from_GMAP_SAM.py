from Bio import Seq
import iCEC
def write_report(err, output_filename):
    nts = ['A','T','C','G']
    with open(output_filename, 'w') as f:
        f.write("--Deletion--\n")
        _sum = sum(err['D'].itervalues())
        for x in nts:
            f.write("{0} -> -: {1:.2f}\n".format(x, err['D'][('-',x)]*1./_sum))
        f.write("--Insertion--\n")
        _sum = sum(err['I'].itervalues())
        for x in nts:
            f.write("{0} -> -: {1:.2f}\n".format(x, err['I'][(x,'-')]*1./_sum))
        f.write("--Substitution--\n")
        _sum = sum(err['S'].itervalues())
        for x in nts:
            for y in nts:
                if x == y: continue
                f.write("{0} -> {1}: {2:.2f}\n".format(x, y, err['S'][(y,x)]*1./_sum))


def main(reader, input, genome, err):
    for r in reader:
        if r.sID not in genome:
            print "skipping", r.sID
            continue
        qstr=Seq.Seq(input[r.qID].sequence.upper())
        tstr=genome[r.sID].seq[r.sStart:r.sEnd].upper()
        if r.flag.strand=='-': qstr=qstr.reverse_complement()
        x=iCEC.unpack_cigar_string(r.cigar)
        iter=iter_cigar(x, qstr, tstr)
        for a,b,c in iter: err[a][(b,c)]+=1

def iter_cigar(cigar_list, qstr, tstr):
    i = 0
    j = 0
    for x in cigar_list:
        if x == 'S':
            i += 1
        elif x == 'N':
            j += 1
        elif x == 'M':
            if qstr[i] != tstr[j]:
                yield 'S',qstr[i], tstr[j]
            i += 1
            j += 1
        elif x == 'D':
            yield 'D','-', tstr[j]
            j += 1
        elif x == 'I':
            yield 'I',qstr[i], '-'
            i += 1
        elif x == 'H' and i > 0:
            return
    print i, j
