# IPython log file

from pbcore.io import BasH5Reader
get_ipython().magic(u'more input.fofn')
reader = BasH5Reader('/home/UNIXHOME/etseng/projects2015/internal_Ting_Maize/runs/BC_1to2kb/0001/Analysis_Results/m141215_033410_42175_c
reader = BasH5Reader('/home/UNIXHOME/etseng/projects2015/internal_Ting_Maize/runs/BC_1to2kb/0001/Analysis_Results/m141215_033410_42175_c100723491910000001823153704301530_s1_p0.1.bax.h5)
reader = BasH5Reader('/home/UNIXHOME/etseng/projects2015/internal_Ting_Maize/runs/BC_1to2kb/0001/Analysis_Results/m141215_033410_42175_c100723491910000001823153704301530_s1_p0.1.bax.h5)
reader = BasH5Reader('/home/UNIXHOME/etseng/projects2015/internal_Ting_Maize/runs/BC_1to2kb/0001/Analysis_Results/m141215_033410_42175_c100723491910000001823153704301530_s1_p0.1.bax.h5')
reader.sequencingZmws
bas = reader[165]
bas
bas.regionTable
9523-9466
reader2 = BasH5Reader('../m141215_033410_42175_c100723491910000001823153704301530_s1_p0.1.ccs.h5')
bas
bas2 = reader2[165]
bas2.numPasses
bas
bas.subreads
bas.adapterRegions
bas.holeNumber
bas.zmwName
len(bas.adapterRegions) - 1
get_ipython().magic(u'ls ')
get_ipython().magic(u'more test_flnc.fasta')
bas = reader[215]
bas2 = reader2[21]
bas2 = reader2[215]
bas.subreads
bas.adapterRegions
len(bas.adapterRegions) - 1
bas2.numPasses
get_ipython().magic(u'ls ')
get_ipython().magic(u'logstart get_num_passes.py')
f = open('num_passes.txt', 'w')
f.write("zmw\tnum_full_pass\tnum_pass\n")
f = open('num_passes.txt', 'w')
f.write("zmw\tnum_full_pass\n')
f.write("zmw\tnum_full_pass\n")
for line in open('input.fofn'):
    reader = BasH5Reader(line.strip())
    for x in reader.sequencingZmws:
        bas = reader[x]
        f.write("{0}\t{1}\n".format(bas.zmw, len(bas.adapterRegions)-1))
        
bas
bas.adapterRegions
bas.zmwName
f = open('num_passes.txt', 'w')
f.write("zmw\tnum_full_pass\n")
for line in open('input.fofn'):
    reader = BasH5Reader(line.strip())
    for x in reader.sequencingZmws:
        bas = reader[x]
        f.write("{0}\t{1}\n".format(bas.zmwName, max(0, len(bas.adapterRegions)-1)))
        
f = open('num_passes.txt', 'w')
f.write("zmw\tnum_full_pass\n")
for line in open('input.fofn'):
    reader = BasH5Reader(line.strip())
    for x in reader.sequencingZmws:
        bas = reader[x]
        n = max(0, len(bas.adapterRegions)-1))
        if n > 0:
            f.write("{0}\t{1}\n".format(bas.zmwName, n))
            
for line in open('input.fofn'):
    reader = BasH5Reader(line.strip())
    for x in reader.sequencingZmws:
        bas = reader[x]
        n = max(0, len(bas.adapterRegions)-1)
        if n > 0:
            f.write("{0}\t{1}\n".format(bas.zmwName, n))
            
f.close()
l
get_ipython().magic(u'ls ')
f.name
get_ipython().magic(u'more num_passes.txt')
get_ipython().magic(u'more num_passes.txt')
