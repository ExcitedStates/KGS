#!/usr/bin/python

import math,sys
import os
import numpy as np
import fnmatch
from pdb_structure import *
from combineKGSPath_steps import extractPath

__description__ = \
"""
Temp factor is adapted according to torsion fluctuations along the path.
Provides pml script that colors according to temp factor.
"""

__author__ = "Dominik Budday"
__date__ = "151011"

import os, sys

def pdbAlterBFactor(pdb,resList,minId):
	# For and ATOM record, update tempFactor number
	out = []
	for line in pdb:
		if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
			# tokens = line.split(' ')
			# print tokens
			resId=int(line[22:26])
			try:
				val=resList[resId-minId][1]
			except IndexError:
				val=0.0
			out.append("%s%6s%s" % (line[0:60],val,line[67:]))
			#print resId, val
		else:
			out.append(line)
	return out

def getTorsionRMSF(pdbPath,pathList,reversePathList,backboneOnly):
	pathFileSepIdx = pdbPath.rfind("/")
	outputDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
	# sampleNum = int(pdbPath[pdbPath.rfind("_")+1:pdbPath.rfind(".pdb")])
	modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
	
	#Identify target name
	included_extensions = ['*[0-9].pdb'] ;
	all_files =[fn for fn in os.listdir(outputDir) if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) ];
	target_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and not fn.startswith(modelName)];
	target_name = target_names[0]
	targetName = str(target_name[target_name.rfind("/")+1:target_name.rfind("_new_")])
	print "targets",len(target_names)
	dofIDs = []
	resids = []
	torsions =[]
	backbone = []

	reverse_dofIDs = []	
	reverse_resids = []
	reverse_torsions =[]
	reverse_backbone = []

	fwdSamples = len(pathList)
	revSamples = len(reversePathList)
	
	# extract information along the forward path
	firstFile = 1
	for sample in pathList:
		with open(outputDir+"/"+modelName+"_q_"+str(sample)+".txt", "r") as q_file:
			if firstFile ==1:
				#Individual information (use first file to extract dof IDs etc.)
				for line in q_file:
					tokens = line.split(' ')
					dofIDs.append(int(tokens[0]))
					resids.append(int(tokens[1]))
					torsions.append(float(tokens[3])*float(tokens[3]))
					backbone.append(float(tokens[5]))
				firstFile=0
			else:
				#Individual torsion information per dof
				count=0
				for line in q_file:
					tokens = line.split(' ')

					#Average
					torsions[count]+=float(tokens[3])*float(tokens[3])

					# #Test the max
					# if float(tokens[3])*float(tokens[3]) > torsions[count]:
					# 	torsions[count]=float(tokens[3])*float(tokens[3])
					count+=1

	# Now we have the squared sum of torsions over the forward path
	#Average over all states along path
	for i in range(len(torsions)):
		torsions[i] = math.sqrt(torsions[i]/fwdSamples)

		# #Test the max (adapt together with top lines)
		# torsions[i] = math.sqrt(torsions[i])

	# extract information along the reverse path
	firstFile = 1
	for sample in reversePathList:
		print sample
		with open(outputDir+"/"+targetName+"_q_"+str(sample)+".txt", "r") as q_file:
			if firstFile ==1:
				#Individual information (use first file to extract dof IDs etc.)
				for line in q_file:
					tokens = line.split(' ')
					reverse_dofIDs.append(int(tokens[0]))
					reverse_resids.append(int(tokens[1]))
					# if(reverse_resids[-1]==99):
					# 	print "DOF value at "+str(reverse_dofIDs[-1])+": "+tokens[3]
					reverse_torsions.append(float(tokens[3])*float(tokens[3]))
					reverse_backbone.append(float(tokens[5]))
				firstFile=0
			else:
				#Individual torsion information per dof
				count=0
				for line in q_file:
					tokens = line.split(' ')

					#Average
					reverse_torsions[count]+=float(tokens[3])*float(tokens[3])

					# #Test the max
					# if float(tokens[3])*float(tokens[3]) > reverse_torsions[count]:
					# 	reverse_torsions[count]=float(tokens[3])*float(tokens[3])

					# if(int(tokens[1])==99):
					# 	print "DOF value at "+tokens[0]+": "+tokens[3]
					count+=1

	# Now we have the squared sum of torsions over the reverse path
	for i in range(len(reverse_torsions)):
		#Average over all states along path
		reverse_torsions[i] = math.sqrt(reverse_torsions[i]/revSamples)

		# #Test the max (adapt together with top lines)
		# reverse_torsions[i] = math.sqrt(reverse_torsions[i])

	return torsions,reverse_torsions
	
	# if desired, also average over the full residue (by number of dofs?)

def main():

	# Procedure
	# extract path ids
	# read-in torsion files
	# compute rmsf on torsion change for each conformation (already given wrt to initial conf)
	# mean for each residue over the forward and reverse path, then average (see rmsf_states)
	# color according to value

	"""
	Outputs pdb and pml file to color according to dihedral rmsf along the path.
	"""

	if len(sys.argv)<4:
		print "Usage: "+sys.argv[0]+" <combinedQ_file>, <backboneOnly? 0:1>, <path.pdb files in a row> "
		sys.exit(1)

	fwdSamples = 0
	revSamples = 0

	resids = []
	dofIDs = []
	torsions =[]
	sumSteps = []
	backbone = []
	
	# Read in residue and backbone information for all dofs
	with open(sys.argv[1]) as outputFile:
		for line in outputFile:
			tokens = line.split(' ')
			dofIDs.append(int(tokens[0]))
			resids.append(int(tokens[1]))
			backbone.append(float(tokens[5]))
			
	backboneOnly=int(sys.argv[2])
	resList=[]
	minId=min(resids)
	maxId=max(resids)
	forwardTorsionRMSF=[0]*len(resids)
	reverseTorsionRMSF=[0]*len(resids)
	overallRMSF = []
	
	for i in range(len(sys.argv)-3): #Go through multiple path files and add them (normally this is only 1 path file)
		pdbPath=sys.argv[i+3]
		pathList, reversePathList = extractPath(pdbPath) #List of path id's
		fwdSamples += len(pathList)
		revSamples += len(reversePathList)
		fwdTorsionRMSF, revTorsionRMSF = getTorsionRMSF(pdbPath, pathList, reversePathList,backboneOnly)
		print len(fwdTorsionRMSF),len(revTorsionRMSF)
		for j in range(len(fwdTorsionRMSF)):
			forwardTorsionRMSF[j] += fwdTorsionRMSF[j]
		for j in range(len(revTorsionRMSF)):
			reverseTorsionRMSF[j] += revTorsionRMSF[j]
	
	numSamples = fwdSamples+revSamples

	#Average on forward and reverse path
	for i in range(len(resids)):
		# Average, forward or reverse
		overallRMSF.append(float((fwdSamples * forwardTorsionRMSF[i] + revSamples * reverseTorsionRMSF[i])/numSamples))
		# overallRMSF.append(reverseTorsionRMSF[i])
		# overallRMSF.append(forwardTorsionRMSF[i])
		# print overallRMSF[i]

	#Prepare list for residue output
	for i in range(minId,maxId+1):
		resVal=[]
		resVal.append(i)
		resVal.append(0)
		resList.append(resVal)

	modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
	
	pdb_file="../"+modelName+".pdb"
	# Read in the pdb file
	f = open(pdb_file,'r')
	pdb = f.readlines()
	f.close()

	# #Now assign maximum DOF values to each residue
	for i in range(len(resids)):
		currentId=resids[i]
		currentVal=overallRMSF[i]

		# if currentVal > 0.2:
		# 	print currentId
		if(backboneOnly==1):
			if(currentVal > resList[currentId-minId][1] and backbone[i]==1): 
				resList[currentId-minId][1]=currentVal
		else:
			if(currentVal > resList[currentId-minId][1]):
				resList[currentId-minId][1]=currentVal
	for i in range(len(resList)):
		val=resList[i][1]
		resList[i][1] = round(val * 100)

	#Mean over residue
	# count = [0]*len(resList)
	# for i in range(len(resids)):
	# 	currentId=resids[i]
	# 	currentVal=overallRMSF[i]
	# 	if(backboneOnly==1):
	# 		if(backbone[i]==1): 
	# 			resList[currentId-minId][1]+=currentVal
	# 			count[currentId-minId]+=1
	# 	else:
	# 		resList[currentId-minId][1]+=currentVal
	# 		count[currentId-minId]+=1

	# for i in range(len(resList)):
	# 	val=resList[i][1]
	# 	resList[i][1] = round(float(val/count[i]) * 100) #to reduce it to non-float ints to keep pdb file layout
		
	# print len(resList)
	out = pdbAlterBFactor(pdb,resList,minId)
	
	d="comboAnalysis"
	if not os.path.exists(d):
		os.makedirs(d)
	os.chdir("comboAnalysis")
	
	fileName = pdb_file[pdb_file.find("/")+1:]
	if backboneOnly==1:
		out_file = "%s_torsionRMSF_bb.pdb" % fileName[:-4]
	else:
		out_file = "%s_torsionRMSF_all.pdb" % fileName[:-4]
	g = open(out_file,'w')
	g.writelines(out)
	g.close()
	
	if backboneOnly==1:
		outPML = "%s_torsionRMSF_bb.pml" % fileName[:-4]
	else: 
		outPML = "%s_torsionRMSF_all.pml" % fileName[:-4]
		
	orig_stdout = sys.stdout
	fout = file(outPML, 'w')
	sys.stdout = fout
	print "from pymol import cmd"
	print "from pymol.cgo import *"

	print "load "+out_file

	print "spectrum b"
	print "bg_color white"
	print "clip slab, 200"
	print "center all"
	print "zoom"

	# hide everything
	print "show lines"
	print "set line_width = 2"
	print "show cartoon"
	# print "set cartoon_transparency, 0.7"
	print "set ray_opaque_background, on"
	print "set cartoon_fancy_helices, on"
	print "cartoon putty"
	print "hide lines"
	# show dashes, hbonds

	sys.stdout = orig_stdout
	fout.close()

if __name__ == "__main__":
    main()
