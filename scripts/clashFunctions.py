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
from math import log
import csv
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import operator

def getClashes(pdbPath,pathList,reversePathList):
	#This function returns forward and reverse clashes separately
	
	# print "Path ini:"+str(pathList[0])+", "+str(pathList[1])

	regPlanner=1
	it_ids=[]
	currentIt=0
	currentPathEntry=0
	pathLength=len(pathList)

	with open("output.txt") as outputFile:
		for line in outputFile:
			if currentPathEntry == pathLength:
				break
			if "Using regular planner" in line:
				regPlanner = 1
			if "Switching to reversePlanner" in line:
				regPlanner = 0
			if "New structure" in line:
				currentIt = currentIt + 1
				
				if regPlanner == 1:
					tokens = line.split(' ')
					num =tokens[3]
					pathId =int(num[num.find('_')+1:num.find('.')])
					if pathId == pathList[currentPathEntry]:
						it_ids.append(currentIt)
						currentPathEntry = currentPathEntry+1
	
	currentClashes = []
	fwdClashes=[]
	atClash=0
	
	currentPathEntry=0
	atConfig=0
	currentIt=1

	# print "All forward confs: "
	# print it_ids
	# print "num path ids "+str(len(pathList))+", num its "+str(len(it_ids))
	
	with open("kgs_planner.log", "r") as log_file:
		for line in log_file:
			if line == '\n':
				currentIt += 1
				continue;
			if(currentIt != it_ids[currentPathEntry]):
				continue;
			if "Using clash constraint for atoms:" in line:
				if atClash==0:
					currentClashes=[]
					atClash=1
					
				tokens = line.split(' ')
				id1=int(tokens[5])
				id2=int(tokens[6])

				currentClashes.append((min(id1,id2),max(id1,id2)))
				
			if "Rejected!" in line:
				atClash=0
				currentClashes=[]
			if "Size of pareto front" in line:
				atClash=0
				currentClashes=[]
			if " .. Overall dofs:" in line:
				if atClash==1:
					# print "Adding clashes at it "+str(currentIt)+": "
					# print currentClashes
					fwdClashes.append(currentClashes)
				
				if currentPathEntry==len(pathList)-1:
					# print "Stopping, as entry "+str(currentPathEntry)+" and "+str(len(pathList)-1)
					break;
				currentPathEntry=currentPathEntry+1
				atClash=0
				# print "Current It "+str(currentIt)+" current path entry "+str(it_ids[currentPathEntry])
				currentClashes=[]
	

	#------------REVERSE-------------Same analysis for the reverse path
	regPlanner=1
	it_ids=[]
	currentIt=0
	currentPathEntry=0
	pathLength=len(reversePathList)

	with open("output.txt") as outputFile:
		for line in outputFile:
			if currentPathEntry == pathLength:
				break
			if "Using regular planner" in line:
				regPlanner = 1
			if "Switching to reversePlanner" in line:
				regPlanner = 0
			if "New structure" in line:
				currentIt = currentIt + 1
				
				if regPlanner == 0:
					tokens = line.split(' ')
					num =tokens[3]
					pathId =int(num[num.find('_')+1:num.find('.')])
					if pathId == reversePathList[currentPathEntry]:
						it_ids.append(currentIt)
						currentPathEntry = currentPathEntry+1
						
	currentClashes = []
	revClashes=[]
	atClash=0
	
	currentPathEntry=0
	atConfig=0
	currentIt=1
	# currentConf=reversePathList[currentPathEntry]

	# print "All reverse confs: "
	# print it_ids

	with open("kgs_planner.log", "r") as log_file:
		for line in log_file:
			if line == '\n':
				currentIt += 1
				continue;
			if(currentIt != it_ids[currentPathEntry]):
				continue;
			if "Using clash constraint for atoms:" in line:
				if atClash==0:
					currentClashes=[]
					atClash=1
					
				tokens = line.split(' ')
				id1=int(tokens[5])
				id2=int(tokens[6])

				currentClashes.append((min(id1,id2),max(id1,id2)))
				
			if "Rejected!" in line:
				atClash=0
				currentClashes=[]
			if "Size of pareto front" in line:
				atClash=0
				currentClashes=[]
			if " .. Overall dofs:" in line:
				if atClash==1:
					# print "Adding clashes at it "+str(currentIt)+": "
					# print currentClashes
					revClashes.append(currentClashes)

				if currentPathEntry==len(reversePathList)-1:
					# print "Stopping, as entry "+str(currentPathEntry)+" and "+str(len(reversePathList)-1)
					break;
				currentPathEntry=currentPathEntry+1
				atClash=0
				# print "Current It "+str(currentIt)+" current path entry "+str(it_ids[currentPathEntry])
				currentClashes=[]

	revClashes.reverse

	return fwdClashes, revClashes;

def getAllClashes(pdbPath,pathList,reversePathList):
	# This function directly combines forward and reverse clashes

	# print "Path ini:"+str(pathList[0])+", "+str(pathList[1])

	regPlanner=1
	it_ids=[]
	currentIt=0
	currentPathEntry=0
	pathLength=len(pathList)

	with open("output.txt") as outputFile:
		for line in outputFile:
			if currentPathEntry == pathLength:
				break
			if "Using regular planner" in line:
				regPlanner = 1
			if "Switching to reversePlanner" in line:
				regPlanner = 0
			if "New structure" in line:
				currentIt = currentIt + 1
				
				if regPlanner == 1:
					tokens = line.split(' ')
					num =tokens[3]
					pathId =int(num[num.find('_')+1:num.find('.')])
					if pathId == pathList[currentPathEntry]:
						it_ids.append(currentIt)
						currentPathEntry = currentPathEntry+1
	
	currentClashes = []
	fwdClashes=[]
	atClash=0
	
	currentPathEntry=0
	atConfig=0
	currentIt=1

	# print "All forward confs: "
	# print it_ids
	# print "num path ids "+str(len(pathList))+", num its "+str(len(it_ids))
	
	with open("kgs_planner.log", "r") as log_file:
		for line in log_file:
			if line == '\n':
				currentIt += 1
				continue;
			if(currentIt != it_ids[currentPathEntry]):
				continue;
			if "Using clash constraint for atoms:" in line:
				if atClash==0:
					currentClashes=[]
					atClash=1
					
				tokens = line.split(' ')
				id1=int(tokens[5])
				id2=int(tokens[6])

				currentClashes.append((min(id1,id2),max(id1,id2)))
				
			if "Rejected!" in line:
				atClash=0
				currentClashes=[]
			if " .. Overall dofs:" in line:
				if atClash==1:
					# print "Adding clashes at it "+str(currentIt)+": "
					# print currentClashes
					fwdClashes.append(currentClashes)

				if currentPathEntry==len(pathList)-1:
					# print "Stopping, as entry "+str(currentPathEntry)+" and "+str(len(pathList)-1)
					break;
				currentPathEntry=currentPathEntry+1
				atClash=0
				# print "Current It "+str(currentIt)+" current path entry "+str(it_ids[currentPathEntry])
				currentClashes=[]
	
	#------------REVERSE-------------Same analysis for the reverse path
	regPlanner=1
	it_ids=[]
	currentIt=0
	currentPathEntry=0
	pathLength=len(reversePathList)

	with open("output.txt") as outputFile:
		for line in outputFile:
			if currentPathEntry == pathLength:
				break
			if "Using regular planner" in line:
				regPlanner = 1
			if "Switching to reversePlanner" in line:
				regPlanner = 0
			if "New structure" in line:
				currentIt = currentIt + 1
				
				if regPlanner == 0:
					tokens = line.split(' ')
					num =tokens[3]
					pathId =int(num[num.find('_')+1:num.find('.')])
					if pathId == reversePathList[currentPathEntry]:
						it_ids.append(currentIt)
						currentPathEntry = currentPathEntry+1
						
	currentClashes = []
	revClashes=[]
	atClash=0
	
	currentPathEntry=0
	atConfig=0
	currentIt=1
	# currentConf=reversePathList[currentPathEntry]

	# print "All reverse confs: "
	# print it_ids

	with open("kgs_planner.log", "r") as log_file:
		for line in log_file:
			if line == '\n':
				currentIt += 1
				continue;
			if(currentIt != it_ids[currentPathEntry]):
				continue;
			if "Using clash constraint for atoms:" in line:
				if atClash==0:
					currentClashes=[]
					atClash=1
					
				tokens = line.split(' ')
				id1=int(tokens[5])
				id2=int(tokens[6])

				currentClashes.append((min(id1,id2),max(id1,id2)))
				
			if "Rejected!" in line:
				atClash=0
				currentClashes=[]
			if " .. Overall dofs:" in line:
				if atClash==1:
					# print "Adding clashes at it "+str(currentIt)+": "
					# print currentClashes
					revClashes.append(currentClashes)

				if currentPathEntry==len(reversePathList)-1:
					# print "Stopping, as entry "+str(currentPathEntry)+" and "+str(len(reversePathList)-1)
					break;
				currentPathEntry=currentPathEntry+1
				atClash=0
				# print "Current It "+str(currentIt)+" current path entry "+str(it_ids[currentPathEntry])
				currentClashes=[]

	revClashes.reverse

	allClashes=fwdClashes
	allClashes.extend(revClashes)

	return allClashes

def getAtomResidueList(pdbFileIn):
	#This function returns a map between atoms and residues for the given pdb file

	atomResidueList={}
	with open(pdbFileIn,"r") as pdbFile:
		for line in pdbFile:
			if "ATOM" == str(line[0:4]):
				atomId=str(line[4:11])
				atomId=int(atomId.lstrip())

				resId=line[22:26]
				resId = int(resId.lstrip())
				atomResidueList[atomId]=resId
	return atomResidueList

def collectAtomClashes(allClashes):
	#This function returns a sorted list of all atom ID-based clashes

	clashCollection={}
	
	for confIt in allClashes:
		for clash in confIt:
			if clash in clashCollection:
				oldVal=clashCollection[clash]
				clashCollection[clash]=oldVal+1
			else:
				clashCollection[clash]=1
				
	sorted_collection = sorted(clashCollection.items(), key=operator.itemgetter(1))
	sorted_collection.reverse()
	return sorted_collection

def collectResidueClashes(clashCollection,atomClashes, atomResidueList):
	# This function returns pairwise residue clashes (not sorted)

	for confIt in atomClashes:
		for clash in confIt:
			atom1=clash[0]
			atom2=clash[1]
			
			resId1 = atomResidueList[atom1]
			resId2 = atomResidueList[atom2]

			resClash = (resId1,resId2)

			if resClash in clashCollection:
				oldVal=clashCollection[resClash]
				clashCollection[resClash]=oldVal+1
			else:
				clashCollection[resClash]=1
				
	# sorted_collection = sorted(clashCollection.items(), key=operator.itemgetter(1))
	# sorted_collection.reverse()
	# return sorted_collection

	return clashCollection

def convertAtomClashesToResidueClashes(clashCollection,atomResidueList,minClashNumber=1,numRuns=1):
	#This function takes a collection of atom-id-based clashes and turns it into pairswise residue clashes
	residueLinks={}

	for entries in clashCollection:
		val = operator.itemgetter(1)(entries)
		key = operator.itemgetter(0)(entries)
		atom1=key[0]
		atom2=key[1]

		val = round(val/numRuns)
		#Atom to residue informaiton
		if val >= minClashNumber: # possibly filter on minimum number of clash occurence (average per run)
			resId1 = atomResidueList[atom1]
			resId2 = atomResidueList[atom2]
			if resId1 == resId2:
				continue
			pairVal=(resId1,resId2)
			if pairVal in residueLinks:
				residueLinks[pairVal] += val
			else:
				residueLinks[pairVal] = val

	return residueLinks

def convertResidueClashesToLinks(clashCollection,minClashNumber=1,numRuns=1):
	# This function filters residue-based clashes by a minimum occurence per run on average 
	residueLinks={}

	for entries in clashCollection:
		val = clashCollection[entries]
		val = round(val/numRuns)
		resId1 = operator.itemgetter(0)(entries)
		resId2 = operator.itemgetter(1)(entries)
		if resId1 == resId2:
			continue
		if val >= minClashNumber:
			residueLinks[entries]=val

	# print residueLinks

	return residueLinks

def convertClashesToResidueNetworks(clashCollection,minClashNumber=1,numRuns=1):
	#This function creates a contiguous network, given a minimum residue-clash occurence per run on average 

	networks=[] #map of contiguous steric networks

	listOfSets=[]

	# Treat minClashNumber as an average over all runs

	for entries in clashCollection:
		val = operator.itemgetter(1)(entries)
		val = round(val/numRuns)
		key = operator.itemgetter(0)(entries)
		resId1=key[0]
		resId2=key[1]
		if resId1 == resId2:
			continue
		#Atom to residue informaiton
		if val >= minClashNumber:
			listOfSets.append(set([resId1,resId2]))

	while True:
		merged_one = False
		supersets = [listOfSets[0]]

		for s in listOfSets[1:]:
			in_super_set = False
			for ss in supersets:
				if s & ss:
					ss |= s
					merged_one = True
					in_super_set = True
					break

			if not in_super_set:
				supersets.append(s)

		#print supersets
		if not merged_one:
			break

		listOfSets = supersets

	supersets.sort(key=len,reverse=True) #sort by length of entries
	residues={} #map of residues and network ID
	numSets = len(supersets)

	id=1
	for net in supersets:
		# print net
		for res in net:
			residues[res]=id
			# print residues[res]
		id=id+1

	# sorted_collection = sorted(residues.items(), key=operator.itemgetter(1))
	# sorted_collection.reverse()

	return residues,numSets

def convertClashesToAllResidues(clashCollection,atomResidueList):
	#This function converts atom-based clashes to residues, accounting for the number of clashes, including residues without a clash
	#This function is mostly used for the pymol output

	residues={}

	for entries in clashCollection:
		val = operator.itemgetter(1)(entries)
		key = operator.itemgetter(0)(entries)
		atom1=key[0]
		atom2=key[1]

		#Atom to residue informaiton
		resId1 = atomResidueList[atom1]
		resId2 = atomResidueList[atom2]

		#Add the first residue
		if resId1 in residues:
			oldVal=residues[resId1]
			residues[resId1]=oldVal+val
		else:
			residues[resId1]=val

		#Add second residue, if not the same as resid1 to prevent double counting of internal clashes
		if resId1 != resId2:
			# if resId2 == 55:
			# 	print "Res 55, atom "+str(atom2)+", clash with res "+str(resId1)+", atom : "+str(atom1)+", "+str(val)+" times."
			if resId2 in residues:
				oldVal=residues[resId2]
				residues[resId2]=oldVal+val
			else:
				residues[resId2]=val

	# Add remaining residues without clashes
	for atom in atomResidueList:
		res = atomResidueList[atom]
		if res in residues:
			continue
		else:
			residues[res]=0

	#Identify maximum value for normalization
	maxVal=0
	for val in residues:
		nums = residues[val]
		if nums > maxVal:
			maxVal=nums

	# Normalize entries
	# for val in residues:
	# 	nums = residues[val]
	# 	residues[val] = round(float(nums)/float(maxVal)*100)/100.0

	sorted_collection = sorted(residues.items(), key=operator.itemgetter(1))
	sorted_collection.reverse()
	return sorted_collection,maxVal

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
			if resId in clashResidues:
				val=clashResidues[resId]
			out.append("%s%6s%s" % (line[0:60],val,line[67:]))
		else:
			out.append(line)
	return out

