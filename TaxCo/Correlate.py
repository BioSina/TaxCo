'''
TaxCo: Correlate

@author: Sina Beier
'''
import re
from scipy import stats as scistats
import numpy
import collections as c
import os

numpy.seterr(all='ignore')

#Correlations for all input levels
def correlate(directory,coefficient,pvalue,cutoff):
    for f in os.listdir(directory):
        if f.endswith(".taxcoIn"):
            correlateFile(f,directory,coefficient,pvalue,cutoff)
 
#Correlation for Taxa           
def correlateFile(f,directory,coefficient,pvalue, cutoff):

    develdict = c.OrderedDict()
    
    infile = open(directory+"/"+f, 'rU')
    line = infile.readline()
    
    for line in infile:
        line = re.sub("\n", "", line)
        split = re.split("\t", line)
        score = split[1:]
        develdict[split[0]] = score
        
    taxondict = c.OrderedDict()
    taxonlist = list()
    
    name1 = re.sub("taxcoIn","taxcoList",f)
    name1 = re.split("\.",name1)
    name1 = name1[0]+"_"+coefficient+"_"+str(pvalue)+"."+name1[1]
    outfile1 = open(directory+"/"+name1,'w')
    outfile1.write("Taxon1\tTaxon2\t"+coefficient+"\tp-value")
    
    usedDict = dict()
        
    for taxon in develdict:
        if not(taxon in taxondict):
            thelist = list()
            taxonlist.append(taxon)
        else:
            thelist = taxondict[taxon]
            
        for taxon2 in develdict:
            if (taxon == taxon2):
                correlation = (0.0, 0.0)
                
            
            else:
                if coefficient== "pearson":
                    scores1 = list(map(float, develdict[taxon]))
                    scores2 = list(map(float, develdict[taxon2]))
                    correlation = scistats.pearsonr(scores1, scores2)
                elif coefficient== "spearman":
                    correlation = scistats.spearmanr(develdict[taxon], develdict[taxon2], 0)
                elif coefficient== "kendall":
                    correlation = scistats.kendalltau(develdict[taxon], develdict[taxon2], True)
                else:
                    print("This correlation coefficient is not supported. Please use spearman, kendall or pearson.")
                    exit(-1)
            if correlation == 1:
                correlation = (numpy.NaN, numpy.NaN)
            rounded = correlation[0]
            rounded = round(abs(rounded),2)
            if((correlation[1]<pvalue) and (rounded > cutoff)):
                if(not ((taxon in usedDict)and(taxon2 in usedDict[taxon]))):
                    outfile1.write("\n"+taxon+"\t"+taxon2+"\t"+str(correlation[0])+"\t"+str(correlation[1]))
                    if(taxon in usedDict):
                        liste = usedDict[taxon]
                    else:
                        liste = list()
                    liste.append(taxon2)
                    usedDict[taxon] = liste 
                    if(taxon2 in usedDict):
                        liste = usedDict[taxon2]
                    else:
                        liste = list()
                    liste.append(taxon)
                    usedDict[taxon2] = liste    
            thelist.append(correlation)
            taxondict[taxon] = thelist
            
    outfile1.close()
    infile.close()
    
    name = re.sub("taxcoIn","taxcoCor",f)
    name = re.split("\.",name)
    name = name[0]+"_"+coefficient+"_"+str(pvalue)+"."+name[1]
    outfile = open(directory+"/"+name,'w')
    line1 = "Taxon"
    for taxon in taxondict:
        line1 = line1+"\t"+taxon
    outfile.write(line1+"\n")
    for taxon in taxondict:
        outfile.write(taxon)
        for correlation in taxondict[taxon]:
            if correlation[1] < pvalue:
                outfile.write("\t"+str(correlation[0]))
            else:
                outfile.write("\t0.0")
        outfile.write("\n")
        
    outfile.close()