import sys
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalNode
from bx.intervals.intersection import IntervalTree
from collections import defaultdict
from csv import DictReader

class GTF:
    def __init__(self, gtf_filename):
        self.gtf_filename = gtf_filename
        self.genome = defaultdict(lambda: IntervalTree()) # chr --> IntervalTree --> (0-start, 1-end, transcript ID)
        self.transcript = defaultdict(lambda: IntervalTree()) # tID --> IntervalTree --> (0-start, 1-end, {'ith': i-th exon, 'eID': exon ID})
        self.exon = defaultdict(lambda: []) # (0start,1end) --> list of (tID, ith-exon, chr)
        
        self.readGTF(self.gtf_filename)
    
    def readGTF(self, filename):
        """
        GTF files
        (0) chr
        (1) annotation source
        (2) type: gene|transcript|CDS|exon|UTR
        (3) 0-based start
        (4) 1-based end
        (5) ignore
        (6) strand: +|-
        (7) phase
        (8) extra stuff (gene ID, transcript ID...) 
        """
        for line in open(filename):
            if line.startswith('#'): continue # header section, ignore
            if len(line.strip()) == 0: continue # some gtf files have blank lines
            raw = line.strip().split('\t')
            chr = raw[0]
            type = raw[2]
            start0, end1 = int(raw[3]), int(raw[4])
            for stuff in raw[8].split('; '):
                _a, _b = stuff.split()
                if _a == "transcript_id": tID = _b[1:-1] # removing quotes ""
                elif _a == "gene_id": gID = _b[1:-1] # removing quotes ""
                
            if type == 'transcript':
                self.genome[chr].insert(start0, end1, tID)
                ith = 0
            elif type == 'exon':
                self.transcript[tID].insert(start0, end1, {'ith':ith,'chr':chr})
                self.exon[(start0,end1)].append((tID, ith, chr))
                ith += 1
        
    def get_exons(self, tID):
        """
        Return the list of intervals for a given tID
        """
        pp = []
        self.transcript[tID].traverse(pp.append)
        return pp

class Coords(GTF):
    def readGTF(self, filename):
        """
        .coords files
        (0) gene name
        (1) chr
        (2) number of exons
        (3) strand
        (4) list of space-separated 0-based start, 1-based end
        """
        for line in open(filename):
            raw = line.strip().split()
            tID = raw[0]
            chr = raw[1]
            ith = 0
            
            if tID in self.transcript:
                print >> sys.stderr, "duplicate tID {0} seen, ignore!".format(tID)
                continue
            
            for i in xrange(4, len(raw), 2):
                start0 = int(raw[i])
                end1 = int(raw[i+1])
                self.genome[chr].insert(start0, end1, tID)
                self.transcript[tID].insert(start0, end1, {'ith':ith, 'chr':chr})
                self.exon[(start0, end1)].append((tID, ith, chr))
                i += 1
        
            
class btabReader:
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self
    
    def next(self):
        return self.read()            
            
    def read(self):
        """
        start-end will be flipped if start > end !!!
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if cur == self.f.tell():
            raise StopIteration, "EOF reached!"
        raw = line.split('\t')
        chr= raw[0]
        seqid = raw[5]
        rStart1 = int(raw[6])
        rEnd1 = int(raw[7])
        i = raw[-2]
        if rStart1 > rEnd1: rStart1, rEnd1 = rEnd1, rStart1
        return {'chr': chr, 'seqid': seqid, 'rStart1': rStart1, 'rEnd1': rEnd1, 'i': i}

class btabBlockReader(btabReader):
    def next(self):
        recs = [self.read()]
        while recs[-1]['i']==recs[0]['i']:
            cur = self.f.tell()
            recs.append(self.read())
        self.f.seek(cur)    
        return recs[:-1]
    
class gmapGFFReader:
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)
        
    def __iter__(self):
        return self
    
    def next(self):
        return self.read()            
            
    def read(self):
        """
        GFF files
        (0) chr
        (1) annotation source
        (2) type: gene|transcript|CDS|exon|UTR
        (3) 0-based start
        (4) 1-based end
        (5) ignore
        (6) strand: +|-
        (7) phase
        (8) extra stuff (gene ID, transcript ID...) 
        
        For gmap output, a series is delimited by '###' line
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration, "EOF reached!!"
        raw = line.strip().split('\t')
        assert raw[2] == 'gene'
        recs = []
        while True:
            line = self.f.readline().strip()
            if line.startswith('###'):
                return recs
            raw = line.split('\t')
            chr = raw[0]
            type = raw[2]
            start1, end1 = int(raw[3]), int(raw[4])
            if type == 'exon':
                for blob in raw[8].split(';'):
                    if blob.startswith('Name='):
                        seqid = blob[5:]
                    elif blob.startswith('Target='):
                        junk, sstart1, ssend1, strand = blob.split()
                        sstart1 = int(sstart1)
                        ssend1 = int(ssend1)                                                               
                recs.append({'chr': chr, 'seqid': seqid, 'rStart1': start1, 'rEnd1': end1, 'sStart1': sstart1, 'sEnd1': ssend1})
            
        
        
              

def btab_reclist_to_interval_list(recs):
    """
    Return chr, list of IntervalNode
    """
    tree = IntervalTree()
    for rec in recs:
        tree.insert(rec['rStart1'], rec['rEnd1'])
    path = []
    tree.traverse(path.append)
    chr = recs[0]['chr']
    return chr, path
    
  
def CompareSimCoordinatesToAlnPath(simCoordinates, alnPath):
    #
    # do silly little dynamic programming to align sets of exons.
    # This could be done in a while loop if there is a 1-1
    # correspondende of exons that overlap, but if multiple overlap,
    # that could cause problems.
    # 
    nAlnExons = len(alnPath)
    nSimExons = len(simCoordinates)
    scoreMat = [[0 for j in xrange(nSimExons+1) ] for i in xrange(nAlnExons+1) ]
    pathMat  = [[0 for j in xrange(nSimExons+1) ] for i in xrange(nAlnExons+1) ]

    diagonal = 0
    up = 1
    left = 2

    for i in xrange(nAlnExons):
        pathMat[i+1][0] = up
    for j in xrange(nSimExons):
        pathMat[0][j+1] = left
    pathMat[0][0] = diagonal
    #return 0

    for i in xrange(nAlnExons):
        for j in xrange(nSimExons):
            overlapScore = 0
            if len(alnPath[i].find(simCoordinates[j].start, simCoordinates[j].end)) > 0: # overlaps!
                overlapScore = 1 # GetOverlapPercent(alnPair, simCoordinates.exonList[j])
                scoreMat[i+1][j+1] = scoreMat[i][j] + overlapScore
                pathMat[i+1][j+1]  = diagonal
            else :
                order = simCoordinates[j].end <= alnPath[i].start #WhichIntervalIsFirst(alnPair, simCoordinates.exonList[j])
                if order:
                    scoreMat[i+1][j+1] = scoreMat[i][j+1] -1 # we penalize ground truth exons that are skipped!!
                    pathMat[i+1][j+1]  = up
                else:             
                    scoreMat[i+1][j+1] = scoreMat[i+1][j] 
                    pathMat[i+1][j+1]  = left

    
    i = nAlnExons    
    j = nSimExons
    matchedExons = []
    while (i > 0 and j > 0):
        if (pathMat[i][j] == diagonal):
            matchedExons.append((i-1,j-1))
            i = i - 1
            j = j - 1
        elif(pathMat[i][j] == left):
            j = j - 1
        else:
            i = i - 1
    matchedExons.reverse()
    return (scoreMat[nAlnExons][nSimExons], matchedExons)  
        
        
def match_transcript(gtf, chr, exon_path):
    """
    exon_tree is an IntervalTree, so it's already sorted
    """
    num_exon = len(exon_path)
    
    #print 'matching transcript for:', exon_path
    
    tIDs = list(set(gtf.genome[chr].find(exon_path[0].start, exon_path[-1].end)))
    best_score, best_matchedExons, best_tID, best_tNum = 0, None, None, None
    for tID in tIDs:
        t_paths = gtf.get_exons(tID)
        score, matchedExons = CompareSimCoordinatesToAlnPath(exon_path, t_paths)
        #print 'matching:', tID, score, matchedExons
        if score > best_score:
            best_tID = tID
            best_tNum = len(t_paths)
            best_score = score
            best_matchedExons = matchedExons
            
    return {'score':best_score, 'matchedExons':best_matchedExons, 'tID': best_tID, 'tID_num_exons': best_tNum}


def categorize_transcript_recovery(info):
    """
    full --- means that every exon in the tID was covered!
    5missX --- means that the assembled one is missing beginning X exons
    3missY --- means that the assembled one is missing ending Y exons
    skipped --- means that the asseembled one is missing some intermediate exons!
    """
    if len(info['matchedExons']) == info['tID_num_exons']: return 'full'
    msg = ''
    if info['matchedExons'][0][0] > 0: msg += '5miss' + str(info['matchedExons'][0][0])
    if info['matchedExons'][-1][0] < info['tID_num_exons']-1: 
        msg += (';' if msg!='' else '') + '3miss' + str(info['tID_num_exons']-1-info['matchedExons'][-1][0])
        
    if msg == '': # must be missing some ground truth exons!
        return 'skipped'
    return msg
    
    

def main(gtf):
    transcript_tally = {}
    for tID in gtf.transcript: 
        transcript_tally[tID] = [0]*len(gtf.get_exons(tID))
    for r in btabBlockReader('sim_gencode_20x_first1000_test2.gmap.tophits.btab'):
        path = btab_reclist_to_interval_list(r)
        info = match_transcript(gtf, r[0]['chr'], path)
        if info['matchedExons'] is None:
            print >> sys.stderr, "Did not find a match for {0}!".format(r[0]['seqid']) 
            continue
        for i, j in info['matchedExons']:
            transcript_tally[info['tID']][i] += 1
    return transcript_tally
    
def main_pasa(gtf):
    pasa_tally = {}
    for tID in gtf.transcript:
        pasa_tally[tID] = [0]*len(gtf.get_exons(tID))
    pasa = GTF('sim_gencode_20x_first1000_test2.pasa_assemblies.denovo_transcript_isoforms.gtf')
    for tID in pasa.transcript:
        path = pasa.get_exons(tID)
        chr = pasa.exon[(path[0].start,path[0].end)][0][2]
        info = match_transcript(gtf, chr, path)
        if info['matchedExons'] is None:
            print >> sys.stderr, "Did not find a match for {0}!".format(tID)
            continue
        for i, j in info['matchedExons']:
            pasa_tally[info['tID']][i] += 1
    return pasa_tally


def eval_gmap(gtf, gmap_filename, input_filename, output_prefix):
    """
    
    Output:
    <output_prefix>.bad --- list of seqids that had no GMAP output or did not match a transcript
    <output_prefix>.report --- 
      <seqid>, <seqlen>, <seqMatchStart>, <seqMatchEnd>, <transcript/gene ID>, <category:full|5missX|3missY|skipped>, <matchedExons>
    """
    from Bio import SeqIO
    fbad = open(output_prefix+'.bad', 'w')
    fbad.write("seqID\tinfo\n")
    fgood = open(output_prefix+'.report', 'w')
    fgood.write("seqID\tseqLen\tseqMatchStart0\tseqMatchEnd1\trefID\tcategory\tmatches\n")
    seqid_missed = [r.id for r in SeqIO.parse(open(input_filename),'fasta')]
    for recs in gmapGFFReader(gmap_filename):
        chr, path = btab_reclist_to_interval_list(recs)
        seqid = recs[0]['seqid']
        seqlen = int(seqid[seqid.rfind('_')+1:]) - int(seqid[seqid.rfind('/')+1:seqid.rfind('_')])
        try:
            seqid_missed.remove(seqid)
        except ValueError: # already removed, ignore?
            pass
        info = match_transcript(gtf, chr, path)
        if info['matchedExons'] is None:
            fbad.write("{0}\tBAD\n".format(seqid))
        else:
            fgood.write("{seqid}\t{seqlen}\t{smstart0}\t{smend1}\t{refID}\t{cat}\t{mat}\n".format(\
                seqid=seqid, seqlen=seqlen, smstart0=min(r['sStart1'] for r in recs)-1, smend1=max(r['sEnd1'] for r in recs),\
                refID=info['tID'], cat=categorize_transcript_recovery(info), mat=info['matchedExons']))
        
    for seqid in seqid_missed:
        fbad.write("{0}\tMISSED\n".format(seqid))
    fbad.close()
    fgood.close()
                                                                      
def make_exon_report(gtf, gmap_report_filename, output_filename):
    """
    Output for each exon:
    <tID>   <exon number 0-based>  <length>  <coverage>   
    """     
    coverage = defaultdict(lambda: defaultdict(lambda: 0)) # tID --> ith-exon --> count           
    for r in DictReader(open(gmap_report_filename), delimiter='\t'):
        tID = r['refID']
        for i,j in eval(r['matches']):
            coverage[tID][i] += 1
    
    f = open(output_filename, 'w')
    for tID in coverage:
        path = gtf.get_exons(tID)
        for ith, exon in enumerate(path):
            f.write("{0}\t{1}\t{2}\t{3}\n".format(tID, ith, exon.end-exon.start, coverage[tID][ith]))
            
    f.close()
        
        
        
        
        
        
        
        
    
    