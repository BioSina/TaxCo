'''
TaxCo: Plot
Needs dependancies: NetworkX and matplotlib

@author: Sina Beier
'''

import re
import os
import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#read in the network-files and run plotting for each
def network(directory):
    filelist = list()
    for f in os.listdir(directory):
        if f.endswith(".taxcoNet"):
            if(not(testFile(directory+"/"+f))):
                filelist.append(f)
                
    print("I have to generate "+str(len(filelist))+" plots, please be patient.")
    graphdir = directory+"/"+"graphs"
    if not os.path.exists(graphdir):
        os.makedirs(graphdir)
    for fi in filelist:
        name = re.sub("\.taxcoNet","",fi)
        levelfile = re.split("__",name)[0]+".taxcoLevel"       
        plotNetwork(directory, graphdir, levelfile,fi)
 
#Test if file has values (not only header)           
def testFile(testfile):
    test = open(testfile,"rU")
    test.readline()
    empty = True
    for line in test:
        empty=False
        break
    test.close()
    return empty
    

#plot the network using Networkx and matplotlib
def plotNetwork(directory, graphdir, levelfile, infile):
    G = nx.Graph()
    G.clear()
    ldict = dict()
    with open(directory+"/"+levelfile, 'rU') as f:
        f.readline() #getting rid of header
        for line in f:
            line = re.sub("\n", "", line)
            split = re.split("\t", line)
            #ldict holds tupel (level,color)
            ldict[split[0]] = (split[1], split[2])

    with open(directory+"/"+infile, 'rU') as f:
        f.readline() #getting rid of header
        for line in f:
            line = re.sub("\n", "", line)
            split = re.split("\t", line)
            n1 = split[0]
            n2 = split[1]
            w = float(split[2])
            if (n1 in ldict.keys()) and (n2 in ldict.keys()):
                G.add_node(n1, size=ldict[n1][0], color=ldict[n1][1])
                G.add_node(n2, size=ldict[n2][0], color=ldict[n2][1])
                G.add_edge(n1, n2, weight=w, color=split[3])
            

    
    #plot the graph
    pos = nx.fruchterman_reingold_layout(G)
    
    for n in G.nodes():
        nlist = [n]
        s = float(G.node[n]['size'])*25
        c = G.node[n]['color']
        l = n
        llist = dict()
        llist[n] = l
        nx.draw_networkx_nodes(G, pos, nlist, node_size=s, node_color=c)
        nx.draw_networkx_labels(G, pos, llist, font_size=10)
    for e in G.edges():
        elist = [e]
        w = float(G.get_edge_data(e[0], e[1])['weight'])*0.5
        c = G.get_edge_data(e[0], e[1])['color']
        nx.draw_networkx_edges(G, pos, elist, width = w, edge_color=c)
        

    plt.axis('off')
    png = re.sub("taxcoNet", "png", infile)
    plt.savefig(graphdir+"/"+png, format='png', dpi=900)
    plt.close()

    #finally also write graph to GraphML
    gml = re.sub("taxcoNet", "GraphML", infile)
    nx.write_graphml(G, graphdir+"/"+gml)
    G.clear()

                
                
                
    