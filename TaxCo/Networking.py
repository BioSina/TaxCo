'''
TaxCo: Networking

@author: Sina Beier
'''
import re
import os

#convert all Correlation lists into the network format
def writeNetworks(directory):
    for f in os.listdir(directory):
        if f.endswith(".taxcoList"):
            prepareNetwork(directory,f)

#produce network file from taxon and metadata correlation files
def prepareNetwork(directory, infile):
    input1 = open(directory+"/"+infile,'rU')
    filename1 = re.sub("taxcoList","taxcoNet",infile)
    filename2 = re.sub("taxcoList","taxcoMeta",infile)
    input2 = open(directory+"/"+filename2,'rU')
    out1 = open(directory+"/"+filename1, 'w')
    
    out1.write("Taxon1\tTaxon2\tweight\tcolor")
    
    collist = ["#F8F8FF", "#BDBDBD", "#848484", "#ADFF2F", "#FFD700", "#FFA500", "#00FFFF"]
    metacolor = 'black'
    levdict = dict()
    metadict = dict()
    line = input1.readline()
    for line in input1:
        line = re.sub("\n","",line)
        if(not(line=="")):
            split = re.split("\t",line)
            tax1 = rename(split[0])
            tax2 = rename(split[1])
            w = weight(split[2])
            c = coloredge(split[2])
            if(not(tax1 in levdict)):
                levdict[tax1] = getLevel(split[0])
            if(not(tax2 in levdict)):
                levdict[tax2] = getLevel(split[1])
        outline = "\n"+tax1+"\t"+tax2+"\t"+str(w)+"\t"+c
        out1.write(outline)
    input1.close()
    line = input2.readline()
    for line in input2:
        line = re.sub("\n","",line)
        if(not(line=="")):
            split = re.split("\t",line)
            tax1 = rename(split[0])
            metavalue = rename(split[1])
            w = weight(split[2])
            c = coloredge(split[2])
            if(not(tax1 in levdict)):
                levdict[tax1] = getLevel(split[0])
            if(not(metavalue in metadict)):
                metadict[metavalue] = metacolor
        outline = "\n"+tax1+"\t"+metavalue+"\t"+str(w)+"\t"+c
        out1.write(outline)
    input2.close()
    out1.close()
    
    filename2 = re.sub(".taxcoList",".taxcoLevel",infile)
    out2 = open(directory+"/"+filename2, 'w')
    out2.write("Taxon\tlevel\tcolor")
    
    for t in levdict:
        out2.write("\n"+t+"\t"+str(levdict[t])+"\t"+collist[levdict[t]])
    for m in metadict:
        out2.write("\n"+m+"\t"+str(1)+"\t"+metacolor)
    out2.close()


def getLevel(name):
    split = re.split(";",name)
    if(split[len(split)-1]=="unclassified"):
        end = False
        levelit = 0
        count=0
        while not(end):
            if (not(split[levelit]=="unclassified")):
                levelit = levelit+1
            else:
                end = True
            if count == len(split)-1:
                end=True
            count = count+1
    else:
        levelit = len(split)
    return levelit

#rename different versions of "unknown"/ "unclassified" taxa to a consistent name
def rename(name):
    split = re.split(";",name)
    newname = ""
    i = 0
    while i < len(split):
        if (not(split[i]=="unclassified")):
            if(not(split[i]=="unknown")):
                newname = split[i]
            else:
                newname = newname+";unknown"
        i = i+1
    if (newname == "")or(newname=="unclassified"):
        newname = "root"
    return newname

#set edge color
def coloredge(correlation):
    col = 'black'
    if(float(correlation)<0.0):
        col='red'
    else:
        col='blue'
    return col

#scale correlation to edge weight
def weight(correlation):
    return abs(float(correlation)*10)