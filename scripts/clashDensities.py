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
from clashFunctions import getAllClashes
from clashFunctions import getAtomResidueList
from clashFunctions import collectAtomClashes
from clashFunctions import convertClashesToAllResidues
from math import log
import csv
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import operator


__description__ = \
"""
Temp factor is adapted according to clash densities recorded during transition.
Script provides pml file to color accordingly.
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

def pdbAlterBFactor(pdb_file,clashResidues):

	# For an ATOM record, update tempFactor number
	f = open(pdb_file,'r')
	pdb = f.readlines()
	f.close()
	out = []
	for line in pdb:
		if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
			# tokens = line.split(' ')
			# print tokens
			resId=int(line[22:26])
			val=0
			for i in range(0,len(clashResidues)):
				res=operator.itemgetter(0)(clashResidues[i])
				if res==resId:
					val=operator.itemgetter(1)(clashResidues[i])
			out.append("%s%6s%s" % (line[0:60],val,line[67:]))
		else:
			out.append(line)
	return out

def main():
	"""
	Adapts color in pdb to clash densities observed along the tree-path of a kgs pdb-file
	"""

	if len(sys.argv)<3:
		print "Usage: "+sys.argv[0]+" <path.pdb files in a row>, <one individual pdb file for atom/residue connection> "
		sys.exit(1)

	samples=[]
	clashConstraints=[]
	rev_clashConstraints=[]
	pdbFile=sys.argv[-1]
	allClashes=[]
	saveDir = os.getcwd()
	
	for i in range(len(sys.argv)-2):
		pdbPath=sys.argv[i+1]

		pathFileSepIdx = pdbPath.find("/output")
		expDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
		pathFileToOpen = pdbPath[pathFileSepIdx+1:] if pathFileSepIdx!=-1 else pdbPath

		os.chdir(expDir)
		#Id's on the configurations on the path, separate for forward and reverse
		pathList, reversePathList = extractPath(pathFileToOpen)
		allClashes.extend( getAllClashes(pathFileToOpen,pathList, reversePathList) )
		os.chdir(saveDir)

	atomResidueList = getAtomResidueList(pdbFile)
	clashCollection = collectAtomClashes(allClashes)
	clashResidues,numSets = convertClashesToAllResidues(clashCollection,atomResidueList)
		
	out = pdbAlterBFactor(pdbFile,clashResidues)
		
	d="comboAnalysis"
	if not os.path.exists(d):
		os.makedirs(d)
	os.chdir("comboAnalysis")
	# d="shortPathAnalysis"
	# if not os.path.exists(d):
	# 	os.makedirs(d)
	# os.chdir("shortPathAnalysis")
	
	fileName = pdbFile[pdbFile.find("/")+1:]
	out_file = "%s_clashDensities.pdb" % fileName[:-4]
	g = open(out_file,'w')
	g.writelines(out)
	g.close()
		
	numClashes=0
	currPercentage=0
	networkPercentage=0

	for i in range(len(clashResidues)):
		numClashes += operator.itemgetter(1)(clashResidues[i])
	print str(numClashes)+" clashes total"


	#Clash percentage boundaries to define different clash levels
	# percentageList = [30,60,80,90,101]
	percentageList = [20,50,70,90,101]
	# percentageList = [10,30,50,80,101]
	colorBorders = []
	percentIndex=0
	for i in range(len(clashResidues)): #range(0,50):
		if (operator.itemgetter(1)(clashResidues[i]) > 0):
			currClashVal = float(operator.itemgetter(1)(clashResidues[i]))
			currPercentage += currClashVal/float(numClashes)
			resPercentage = float(i+1)/float(len(clashResidues))
			
			print str(clashResidues[i])+" "+str(round(resPercentage*10000)/100)+" percent res, "+str(round(currPercentage*10000)/100)+" percent clashes"
			if round(currPercentage*10000)/100 >= percentageList[percentIndex]:
				colorBorders.append(currClashVal-0.5)
				percentIndex += 1

	sys.stdout.write("select clashingResidues, resi ")
	for i in range(len(clashResidues)):
		if operator.itemgetter(1)(clashResidues[i]) > 0:
			res=operator.itemgetter(0)(clashResidues[i])
			sys.stdout.write(str(res)+"+")
		else: 
			break
	print "Done"

	#Provide pml file for adapted coloring
	outPML = "%s_clashDensities_"+str(percentageList[0])+str(percentageList[1])+str(percentageList[2])+str(percentageList[3])+".pml" % fileName[:-4]
			
	orig_stdout = sys.stdout
	fout = file(outPML, 'w')
	sys.stdout = fout
	print "from pymol import cmd"
	print "from pymol.cgo import *"

	print "load "+out_file

	#Select colors according to clashes
	print "color gray80, b<0.99"
	print "color blue, ( b > 0.99 and b < "+str(colorBorders[-1])+")"
	print "color deepteal, ( b > "+str(colorBorders[-1])+" and b < "+str(colorBorders[-2])+")"
	print "color green, ( b > "+str(colorBorders[-2])+" and b < "+str(colorBorders[-3])+")"
	print "color salmon, ( b > "+str(colorBorders[-3])+" and b < "+str(colorBorders[-4])+")"
	print "color red, ( b > "+str(colorBorders[-4])+" )"
	print "select high, b>"+str(colorBorders[-2])
	print "show sticks, high"

	print "bg_color white"

	# hide everything
	print "show lines"
	print "set line_width = 1"
	print "hide lines, b<0.99"
	print "show cartoon"
	print "set ray_opaque_background, on"
	print "set cartoon_fancy_helices, on"

	sys.stdout = orig_stdout
	fout.close()

if __name__ == "__main__":
	main()
