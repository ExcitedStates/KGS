#!/usr/bin/python
import sys
import os
import numpy as np
import fnmatch
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
mpl.use('Agg')
from combineKGSPath_steps import extractPath
from clashFunctions import getClashes
from clashFunctions import getAtomResidueList
from clashFunctions import collectAtomClashes
from clashFunctions import collectResidueClashes
from clashFunctions import convertAtomClashesToResidueClashes
from clashFunctions import convertClashesToResidueNetworks
from clashFunctions import pdbAlterBFactor
from clashFunctions import convertResidueClashesToLinks
from math import log
import csv
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import operator
import networkx as nx, re


__description__ = \
"""
This script identifies contiguous networks of steric clashes to identify sterically coupled regions in
the molecule
"""

__author__ = "Dominik Budday"
__date__ = "160301"


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# from pylab import rcParams
# rcParams['figure.figsize'] = 6.5/2, 6.5*3/8

ftSize=19
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

def build_graph(paths, network_length_threshold = 0):
    """
    input is a list of lists from extract_pathway
    will connect networks from that list and throw out any networks with less than network_length_threshold nodes
    output is a graph with subgraphs representing different networks, edges are weighted by number of pathways connecting nodes
    """
    pathway_lengths = []
    G=nx.Graph()
    maxVal = 0
    for connection in paths:
    	# print connection
    	weightVal = paths[connection]
    	pathway_exists = 0
    	for e in G.edges():
    	    if connection[0] in e and connection[1] in e:
    	        G[connection[0]][connection[1]]['weight'] += weightVal
    	        if G[connection[0]][connection[1]]['weight'] > maxVal:
    	        	maxVal = G[connection[0]][connection[1]]['weight']
    	        pathway_exists = 1
    	if not pathway_exists:
    	    G.add_edge(connection[0],connection[1],weight=weightVal)
    	    if G[connection[0]][connection[1]]['weight'] > maxVal:
    	    	maxVal = G[connection[0]][connection[1]]['weight']

    connects = nx.connected_components(G)
    ###FILTER GRAPH TO ONLY NETWORKS > THRESHOLD
    # for network in connects:
    #     if len(network) < network_length_threshold:
    #         G.remove_node(network[0])
    #         G.remove_node(network[1])
    
    #for e in G.edges():
#         print e, G.get_edge_data(e[0],e[1])
#     plt.hist(pathway_lengths)
#     plt.show()
#     print len(pathway_lengths) # total number of pathways
    return G, maxVal


def main():
    """
    Adapts b-factor to color according to steric clash networks identified along the tree-path of a kgs pdb-file
    """
    
    if len(sys.argv)<4:
        print "Usage: "+sys.argv[0]+"<minClashNumber>, <output.txt> <path.pdb files in a row> "#, <reverse pdb file>, <forward pdb file> "
        print "Start this from the base directory of all experiments"
        sys.exit(1)
    
    minClashNumber=int(sys.argv[1])
        
    pdbFile = ""
    pdbFileRev = ""
    outputTxtFile = sys.argv[2]
    with open(outputTxtFile) as outputFile:
        for line in outputFile:
            if "--initial " in line:
                pdbFile = line[line.find("--init")+10:line.rfind(".pdb")+4]
                modelName = line[line.rfind("/")+1:line.rfind(".pdb")]
            if "--target " in line:
                pdbFileRev = line[line.find("--target")+9:line.rfind(".pdb")+4]
                break;
            

    # outputPDBDir = outputTxtFile[0:outputTxtFile.rfind("/")]

	# pdbPath=sys.argv[3]
	# if( len(sys.argv) > 4):
	# 	pdbFile=sys.argv[-1]
	# 	modelName = str(pdbFile[pdbFile.rfind("/")+1:pdbFile.rfind(".pdb")])
	# 	pdbFileRev = sys.argv[-2]
	# else:
	# 	modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
	# 	pdbFile = "../"+modelName+".pdb"
	# 	
	# print modelName


    fwdClashes = []
    revClashes = []
    
    currDir = os.getcwd()
    sumRuns=0
    
    # Removed multi-path pdb file support
    # for pFile in range(len(sys.argv)-4):
    # pdbPath=sys.argv[pFile+2]
    for pdbPath in sys.argv[3:]:
        pathFileSepIdx = pdbPath.find("/output")
        expDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
        pathFileToOpen = pdbPath[pathFileSepIdx+1:] if pathFileSepIdx!=-1 else pdbPath
    
        os.chdir(expDir)
        #Id's on the configurations on the path, separate for forward and reverse
        pathList, reversePathList = extractPath(pathFileToOpen)
        # allClashes.extend( getAllClashes(pathFileToOpen,pathList, reversePathList) )
        forwardClashes, reverseClashes = getClashes(pathFileToOpen,pathList, reversePathList)
        fwdClashes.extend(forwardClashes)
        revClashes.extend(reverseClashes)
        sumRuns +=1
        os.chdir(currDir)
    
    #END of multi-path loop
    # os.chdir(outputPDBDir)
    fwdAtomResidueList = getAtomResidueList(pdbFile)
    revAtomResidueList = getAtomResidueList(pdbFileRev)

    
    # This is on an atom-clash based level
    # clashCollection = collectAtomClashes(allClashes)
    # residueLinks = convertAtomClashesToResidueClashes(clashCollection,atomResidueList,minClashNumber)
    
    # This is on an residue-clash based level
    clashCollection = {}
    clashCollection = collectResidueClashes(clashCollection,fwdClashes,fwdAtomResidueList)
    clashCollection = collectResidueClashes(clashCollection,revClashes,revAtomResidueList)
    sorted_collection = sorted(clashCollection.items(), key=operator.itemgetter(1))
    sorted_collection.reverse()
    
    residueLinks = convertResidueClashesToLinks(clashCollection,minClashNumber,sumRuns)
    print residueLinks
    G, maxVal = build_graph(residueLinks)

    #make pictures and output networks
    colorVals = ['red','cyan','magenta','yellow','green','orange','purple','blue','black','gray','forest']
    #Use colorData to specify special residues that are colored differently (e.g. show up in experimental data)
    colorData={}
    # colorData = {55: 1584, 57:1511, 61:1354, 63:1957, 92:1072, 99:1320, 113:1400, 122:1408, 126:1453, 128:1038} #62:1542,
    minColVal = 500
    maxColVal = 1000
    jet = cm = plt.get_cmap('cool')
    cNorm  = colors.Normalize(vmin=minColVal, vmax=maxColVal)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    i = 0
    colorMat={}
    
    nodeSize = 1000
    
    # print maxVal, nodeSize
    
    os.chdir(currDir)
    d="comboAnalysis"
    if not os.path.exists(d):
        os.makedirs(d)
    os.chdir("comboAnalysis")
    
    rounds=0
    C=nx.connected_component_subgraphs(G)
    for g in C:
        resis = []
        for resi in g:
            resis.append(str(resi))
            if resi in colorData:
                # print "Found resi "+str(resi)
                colorMat[resi]=scalarMap.to_rgba(colorData[resi])
            # 	netResis.append(resi)
            # else:
            # 	colorMat[resi] = [160,160,160]
    
        if len(resis) > 2:
            pos=nx.spring_layout(g,weight=0.1,k=0.15,iterations=500,scale=0.1)
            #pos=nx.spring_layout(g,weight=0.1,k=0.4,iterations=1000)
            for res in g:
                if res in colorData:
                    nx.draw_networkx_nodes(g,pos,nodelist=[res],node_size=nodeSize,node_color='blue')
                else:
                    nx.draw_networkx_nodes(g,pos,nodelist=[res],node_size=nodeSize,node_color='grey')
            # nx.draw_networkx_nodes(g,pos,node_size=nodeSize,node_color=colorVals[i])
    
            nx.draw_networkx_labels(g,pos,font_size=ftSize,font_family='sans-serif',font='Arial',font_color='white')
            for e in g.edges():
    #             print e, G.get_edge_data(e[0],e[1])['weight']
                nx.draw_networkx_edges(G,pos,edgelist=[e],width=round(float(G.get_edge_data(e[0],e[1])['weight'])/maxVal/2/minClashNumber*(nodeSize-1)+1))
            print "color %s, resi %s" %(colorVals[i],"+".join(resis))
            print "create %s_sec, resi %s" %(colorVals[i],"+".join(resis))
            print "show surface, %s_sec" %(colorVals[i])
            plt.savefig("%s_%s.png" %(colorVals[i],rounds),dpi=600)
            plt.clf()    
            i+=1
            if i > 8:
                i = 0
                rounds += 1

if __name__ == "__main__":
	main()