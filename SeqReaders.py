from pbcore.io.FastaIO import FastaRecord

class LazyFastaReader:
    """
    NOTE: this version works with pbcore, not Biopython
    
    This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
    are too big to fit entirely to memory. The only requirement is that every id line
    begins with the symbol >. It is ok for the sequences to stretch multiple lines.
    The sequences, when read, are returned as FastaRecord objects.

    Example:
        r = FastaReader('output/test.fna')
        r['6C_49273_NC_008578/2259031-2259297'] ==> this shows the FastaRecord
    """
    def __init__(self, fasta_filename):
        self.f = open(fasta_filename)
        self.d = {}
        
        while 1:
            line = self.f.readline()
            if len(line) == 0: break
            if line.startswith('>'):
                id = line.strip()[1:].split(None, 1)[0] # the header MUST be just 1 line
                if id in self.d:
                    raise Exception, "Duplicate id {0}!!".format(id)
                self.d[id] = self.f.tell()

    def __getitem__(self, k):
        if k not in self.d:
            raise Exception, "key {0} not in dictionary!".format(k)
        self.f.seek(self.d[k])
        content = ''
        for line in self.f:
            if line.startswith('>'):
                break
            content += line.strip()
        return FastaRecord(k, content)

    def keys(self):
        return self.d.keys()
