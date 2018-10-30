'''
TaxCo: Convert

@author: Sina Beier

This is the conversion part of TaxCo, converting input files from QIIME, mothur or MEGAN into the 
"taxcoIn" file format.
'''
import re
import os
import collections as c
import sys

leveled = c.OrderedDict()
corrected = c.OrderedDict()

#Conversion for QIIME
def QIIMEConvert(L2, L3, L4, L5 ,L6 , outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    convertLevel(L2, outdir+"/L2.taxcoIn")
    convertLevel(L3, outdir+"/L3.taxcoIn")
    convertLevel(L4, outdir+"/L4.taxcoIn")
    convertLevel(L5, outdir+"/L5.taxcoIn")
    convertLevel(L6, outdir+"/L6.taxcoIn")
    
#Conversion for QIIME per Level  
def convertLevel(inputfile, outputfile):
    infile = open(inputfile, 'rU')
    outfile = open(outputfile, 'w')
    line = infile.readline()
    outfile.write(line)
    for line in infile:
        split = re.split("\t", line)
        values = split[1:]
        values = "\t".join(values)
        pre_tax = split[0]
        split_tax = re.split(";", pre_tax)
        taxstring = ""
        for i in split_tax:
            name = ""
            if (re.search("__", i)):
                splitit = re.split("__", i)
                if (not(splitit[1]=="")):
                    name = splitit[1]
                else:
                    name = "unknown"
            elif (re.match("Other", i)):
                name = "unclassified"
            name = re.sub("\s","_", name)
            taxstring = taxstring+";"+name
        taxstring = re.sub('^[;]+',"",taxstring)
        outfile.write(taxstring+"\t"+values)
    infile.close()
    outfile.close()

#Conversion for mothur
def mothurConvert(mothur, outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    try:    
        infile = open(mothur,'rU')
    except:
        sys.stderr.write("Input file could not be opened.")
        exit(1)
        
    header = infile.readline()
    header = re.sub("\n","", header)
    header = re.split("\t",header)
    newhead = "Taxon\t"+"\t".join(header[5:])
    for line in infile:
        line = re.sub("\n","",line)
        if(not(line == "")):
            split = re.split("\t",line)
            values = split[5:]
            if(values[len(values)-1]==""):
                values.pop(len(values)-1)
            level = split[0]
            if (not level in leveled):
                leveled[level] = c.OrderedDict()
                corrected[level] = c.OrderedDict()
    
            ranking = split[1]
            name = split[2]
            total = split[4]
            if(not(ranking == "")):
                leveled[level][ranking] = [name,total,values]
    infile.close()
    
    for i in leveled:
        if (not i=="6"):
            for r in leveled[i]:
                suby = subtotal(i,r)
                corrected[i][r] = [leveled[i][r][0],suby[0],suby[1]]
        corrected["6"] = leveled["6"] 

    writeLevel(newhead, leveled, corrected, 2, outdir+"/L2.taxcoIn")
    writeLevel(newhead, leveled, corrected, 3, outdir+"/L3.taxcoIn")
    writeLevel(newhead, leveled, corrected, 4, outdir+"/L4.taxcoIn")
    writeLevel(newhead, leveled, corrected, 5, outdir+"/L5.taxcoIn")
    writeLevel(newhead, leveled, corrected, 6, outdir+"/L6.taxcoIn")
                         
#Subtotal of abundance assigned to lower levels
def subtotal(level, rankID):
    total = leveled[level][rankID][1]
    sub = []
    for i in leveled[level][rankID][2]:
        sub.append(i)
    sublevel = leveled[str(int(level)+1)]
    for rID in sublevel:
        if(re.match('^'+rankID+'\.',rID)):
            total = str(int(total) - int(sublevel[rID][1]))
            vals = sublevel[rID][2]
            j = 0
            while (j<len(vals)):
                sub[j] = str(int(sub[j]) - int(vals[j]))
                j = j+1
    return (total,sub)

def leveledname(level, leveled, rankID):
    full = fullname(leveled,rankID)
    split = re.split(";",full)
    part = split[:level]
    name = ";".join(part)
    return name

def fullname(leveled, rankID):
    name = ""
    ending = ""
    split=re.split("\.",rankID)
    i = 7-len(split)
    while(i > 0):
        ending = ending+";Unclassified"
        i = i-1
    j = 1
    while (j <= len(split)):
        subid = ".".join(split[:j])
        if(not(subid=="0")):
            part = leveled[str(j-1)][str(subid)][0]
            name = name+";"+part
        j = j+1
    name = name+ending
    name = re.sub("^;","",name)
    return name
    
#mothur conversion per level
def writeLevel(header, leveled, corrected, level, outfile):
    output = open(outfile, 'w')
    output.write(header+"\n")
    
    l = 0
    while (l < level):
        strl = str(l)
        for r in corrected[strl]:
            values = "\t".join(corrected[strl][r][2])
            if(not(corrected[strl][r][1])=="0"):
                output.write(leveledname(l, leveled,r)+"\t"+values+"\n")
        l = l+1  
          
    for r in leveled[str(level)]:
        vals = "\t".join(leveled[str(level)][r][2])
        output.write(leveledname(level,leveled,r)+"\t"+vals+"\n")
        
 
#MEGAN conversion   
def MEGANConvert(L2, L3, L4, L5 ,L6 , outdir):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    levelDict  = dict()
    header, levelDict[2] = readMEGAN(L2)
    header, levelDict[3] = readMEGAN(L3)
    header, levelDict[4] = readMEGAN(L4)
    header, levelDict[5] = readMEGAN(L5)
    header, levelDict[6] = readMEGAN(L6)
    
    
    for x in levelDict:
        outfile = open(outdir+"/L"+str(x)+".taxcoIn", 'w')
        outfile.write(header)
        for i in levelDict[x]:
            v = "\t".join(levelDict[x][i])
            outfile.write(i+"\t"+v)
        outfile.close()
    
    
def readMEGAN(infile):
    infile = open(infile, 'rU')
    outdict = c.OrderedDict()
    header = infile.readline()

    for line in infile:
        split = re.split("\t", line)
        values = split[1:]
        pre_tax = split[0]
        pre_tax = re.sub('\"', "", pre_tax)
        split_tax = re.split(";", pre_tax)
        taxstring = split_tax[-2]
        name = re.sub("\s","_", taxstring)
        outdict[name]=values
    infile.close()
    return header, outdict
    