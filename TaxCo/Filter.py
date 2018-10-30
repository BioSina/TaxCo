'''
TaxCo: Filter

@author: Sina Beier
'''

import re
import os
import collections

#filter by given cutoffs for correlation
def filterAll(directory,cutoffs):
    files = list()
    for f in os.listdir(directory):
        if f.endswith(".taxcoNet"):
            files.append(directory+"/"+f)
    for fi in files:
        if isinstance(cutoffs, collections.Iterable):
            for c in cutoffs:
                filterCutoff(fi, float(c))
        else:
            filterCutoff(fi, float(cutoffs))
    for f2 in os.listdir(directory):
        if f2.endswith(".taxcoNet"):
            pos_neg(directory+"/"+f2)
        
    
#write out separate positive and negative correlation files (for plotting)
def pos_neg(infile):
    input1 = open(infile,'rU')
    posname = re.sub(".taxcoNet","__positive.taxcoNet",infile)
    negname = re.sub(".taxcoNet","__negative.taxcoNet",infile)
    pos = open(posname, 'w')
    neg = open(negname, 'w')
    header = input1.readline()
    pos.write(header)
    neg.write(header)
    for line in input1:
        newline = re.sub("\n","",line)
        split = re.split("\t",newline)
        if split[3]=='red':
            neg.write(line)
        if split[3]=='blue':
            pos.write(line)
    input1.close()
    pos.close()
    neg.close()
  
#filter out correlations by cutoff for absolute correlation  
def filterCutoff(infile, cutoff):
    input1 = open(infile,'rU')
    name = re.sub(".taxcoNet","__"+str(cutoff)+".taxcoNet",infile)
    outfile = open(name, 'w')
    header = input1.readline()
    header = re.sub("\n", "", header)
    outfile.write(header)
    for line in input1:
        newline = re.sub("\n", "", line)
        split = re.split("\t",newline)
        if round(abs(float(split[2])),2)>(cutoff*10):
            outfile.write("\n"+newline)
    input1.close()
    outfile.close()