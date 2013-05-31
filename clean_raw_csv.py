import fnmatch
import os
from csv import DictReader

files=fnmatch.filter(os.listdir('.'),'filtered_summary.csv')        
for file in files:
    f=open(file+'.cleaned.csv','w')
    for r in DictReader(open(file),delimiter=','):
        if r['PassedFilter']=='1': f.write(r['Readlength']+'\n')
        
f.close()


files=fnmatch.filter(os.listdir('.'),'*subread*.csv')        
for file in files:
    f=open(file+'.cleaned.csv','w')
    for r in DictReader(open(file),delimiter=','):
        if r['PassedFilter']=='1': f.write(r['Length']+'\n')
        
f.close()




