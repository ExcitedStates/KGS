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
from clashFunctions import getAllClashes
from clashFunctions import getAtomResidueList
from clashFunctions import collectResidueClashes
from clashFunctions import convertClashesToResidueNetworks
from clashFunctions import pdbAlterBFactor
from math import log
import csv
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import operator


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

from pylab import rcParams
rcParams['figure.figsize'] = 6.5/2, 6.5*3/8

ftSize=12
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

def main():

	"""
	Adapts b-factor to color according to steric clash networks identified along the tree-path of a kgs pdb-file
	"""

	if len(sys.argv)<5:
		print "Usage: "+sys.argv[0]+"<minClashNumber>, <path.pdb files in a row>, <reverse pdb file>, <forward pdb file> "
		print "Start this from the base directory of all experiments"
		sys.exit(1)
	
	pdbFile = ""
	pdbFileRev = ""
	pdbPath=sys.argv[2]
	if( len(sys.argv) > 4):
		pdbFile=sys.argv[-1]
		modelName = str(pdbFile[pdbFile.rfind("/")+1:pdbFile.rfind(".pdb")])
		pdbFileRev = sys.argv[-2]
	else:
		modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
		pdbFile = "../"+modelName+".pdb"
		
	# pdbPath = sys.argv[2]
	# allClashes = []
	fwdClashes = []
	revClashes = []
	minClashNumber=int(sys.argv[1])

	currDir = os.getcwd()
	sumRuns=0
	
	# Removed multi-path pdb file support
	for pFile in range(len(sys.argv)-4):
		pdbPath=sys.argv[pFile+2]
	
		pathFileSepIdx = pdbPath.find("/output")
		expDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
		pathFileToOpen = pdbPath[pathFileSepIdx+1:] if pathFileSepIdx!=-1 else pdbPath
	
		print "Changing to "+str(expDir)
		os.chdir(expDir)
		#Id's on the configurations on the path, separate for forward and reverse
		# pathList, reversePathList = extractPath(pdbPath)
		# allClashes = getClashes(pdbPath,pathList, reversePathList)
	
		print "Now in "+str(os.getcwd())
		pathList, reversePathList = extractPath(pathFileToOpen)
		# allClashes.extend( getAllClashes(pathFileToOpen,pathList, reversePathList) )
		forwardClashes, reverseClashes = getClashes(pathFileToOpen,pathList, reversePathList)
		fwdClashes.extend(forwardClashes)
		revClashes.extend(reverseClashes)
		sumRuns += 1
		os.chdir(currDir)
	# End of multi-path loop

	fwdAtomResidueList = getAtomResidueList(pdbFile)
	revAtomResidueList = getAtomResidueList(pdbFileRev)

	# This is on a residue-clash based level
	clashCollection = {}
	clashCollection = collectResidueClashes(clashCollection,fwdClashes,fwdAtomResidueList)
	clashCollection = collectResidueClashes(clashCollection,revClashes,revAtomResidueList)
	sorted_collection = sorted(clashCollection.items(), key=operator.itemgetter(1))
	sorted_collection.reverse()

	clashResidues,numSets = convertClashesToResidueNetworks(sorted_collection,minClashNumber,sumRuns)

	out = pdbAlterBFactor(pdbFile,clashResidues)
	
	print "Number of sets: "+str(numSets)
	d="comboAnalysis"
	if not os.path.exists(d):
		os.makedirs(d)
	os.chdir("comboAnalysis")

	out_file = "%s_clashNetworks_%s.pdb" %(modelName,str(minClashNumber))
	g = open(out_file,'w')
	g.writelines(out)
	g.close()

	jet = cm = plt.get_cmap('jet')
	cNorm  = colors.Normalize(vmin=1, vmax=numSets)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

	#Provide pml file for adapted coloring
	outPML = "%s_clashNetworks_%s.pml" %(modelName,str(minClashNumber))
		
	orig_stdout = sys.stdout
	fout = file(outPML, 'w')
	sys.stdout = fout
	print "from pymol import cmd"
	print "from pymol.cgo import *"

	print "load "+out_file
	print "bg_color white"

	# print "spectrum b"
	# print "color gray80, b<0.99"
	# for i in range(1,numSets):
	# 	colorval=scalarMap.to_rgba(i)
	# 	print "color ["+str(colorval[0])+","+str(colorval[1])+","+str(colorval[2])+"], (b > "+str(float(i)-0.01)+" and b < "+str(float(i)+0.01)+")"

	# List of available colours
	colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
       'salmon', 'smudge', 'hotpink', 'slate', 'magenta', 'orange', 'marine', \
       'olive', 'density', 'teal', 'forest', 'firebrick', 'chocolate',    \
       'wheat', 'white']
	ncolours = len(colours)

	# print "create NW0,b<0.99"
	print "color gray80, b<0.99"
	j=0
	for i in range(1,numSets+1):
		colorval=scalarMap.to_rgba(i)
		print "create NW"+str(i)+", (b > "+str(float(i)-0.01)+" and b < "+str(float(i)+0.01)+")"
		print "hide everything, NW"+str(i)
		print "color "+str(colours[j])+",(b > "+str(float(i)-0.01)+" and b < "+str(float(i)+0.01)+")"
		print "show surface, NW"+str(i)
		j=j+1;
		if(j==ncolours):
			j=0

	# hide everything
	print "show lines"
	print "hide lines, b<0.99"
	print "set line_width = 1"
	print "show cartoon"
	print "set ray_opaque_background, on"
	print "set cartoon_fancy_helices, on"
	print "set transparency, 0.7"

	sys.stdout = orig_stdout
	fout.close()

	os.chdir(currDir)

if __name__ == "__main__":
 	main()