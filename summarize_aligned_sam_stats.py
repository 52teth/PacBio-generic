import os, sys

def summarize(dirname):
    os.chdir(dirname)
    
    refhits = {}
    num_subreads = int(os.popen("grep -c \">\" data_fromMan/filtered_subreads.fasta").read())
    num_aligned = 0
    with open("data/aligned_reads.sam") as f:
        for line in f:
            if line.startswith('@SQ'):
                rid = line.split('\t')[1][3:]
                refhits[rid] = 0
            elif not line.startswith('@'):
                num_aligned += 1
                rid = line.split('\t')[2]
                refhits[rid] += 1
    
    print("Number of subreads: {0}".format(num_subreads))
    print("Number of subreads aligned: {0}".format(num_aligned))
    print("###")
    for rid,hits in refhits.iteritems():
        print("{0}\t{1}".format(rid, hits))
        
        
if __name__ == "__main__":
    summarize(sys.argv[1])
                
                
                
            
    
    
    
    