#!/usr/bin/python
import sys
import numpy as py
from math import floor
import fnmatch
import os

def extractPath(pdbPath):
	reversePathList=[]
	print "Opening "+str(pdbPath)
	fp = open(pdbPath)
	fwd=0
	rev=0
	i=0
	for line in fp:
		if "REMARK	Path is" in line:
			print "Found forward path"
			pathStr = line[line.find("is")+2:].strip()
		if "REMARK	Tree-path" in line:
			pathStr = line[line.find("=")+2:].strip()
			fwd=1
		if "Reverse path follows" in line:
			rev=1
			print "Found reverse path"
			reversePathStr = line[line.find("follows")+8:].strip()
			break;
	
	fp.close()
	if fwd==0:
		pathList = pathStr.split(", ")
	else:
		pathList = pathStr.split(" ")
	if rev==1:
		reversePathList = reversePathStr.split(", ")
		
	del pathList[-1] #Remove last element
	pathList = map(int, pathList) #Convert to list-of-int
	pathList = sorted(pathList) #First sample is first
	#print pathList
	
	if rev==1:
		del reversePathList[-1]
		reversePathList = map(int, reversePathList) #Convert to list-of-int
		reversePathList = sorted(reversePathList) #Path follows from reverse initial
		
	#We return forward and reverse path list, each in order of their individual sampling
	#To combine paths, the reverse path needs to be "reversed"
	print "Forward is "+str(len(pathList))+" samples, reverse is "+str(len(reversePathList))+" samples, overall is "+str(len(pathList)+len(reversePathList))+" samples"
	# print pathList
	# print reversePathList
	return pathList,reversePathList


def writePath(pdbPaths, stepNum=None):
	if stepNum==None:
		stepNum = 10
	else:
		stepNum = int(stepNum)
	
	included_extensions = ['*[0-9].pdb'] ;

	i=1
	ir=1
	for pdbPath in pdbPaths:
		pathList, reversePathList = extractPath(pdbPath)
	
		pathFileSepIdx = pdbPath.rfind("/")
		outputDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
		# sampleNum = int(pdbPath[pdbPath.rfind("_")+1:pdbPath.rfind(".pdb")])
		modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
			
		all_files =[fn for fn in os.listdir(outputDir) if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) ];
		target_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and not fn.startswith(modelName)];
		target_name = target_names[0]
		targetName = str(target_name[target_name.rfind("/")+1:target_name.rfind("_new_")])

		# pathSeq = pathList[::stepNum]
		# reversePathSeq = reversePathList[::stepNum]
		indices=py.linspace(0,len(pathList)-1,num=200)
		for i in range(len(indices)):
			indices[i]=floor(indices[i])
	
		pathSeq = pathList[::stepNum]
		
		indices=py.linspace(0,len(reversePathList)-1)
		for i in range(len(indices)):
			indices[i]=floor(indices[i])
			
		reversePathSeq = reversePathList[::stepNum]

		# Change to comboAnalsis
		currentPath=os.getcwd()
		d="comboAnalysis"
		if not os.path.exists(d):
		    os.makedirs(d)
		os.chdir(d)

		orig_stdout = sys.stdout
		fout = file(modelName+'_path_'+str(stepNum)+'.pdb', 'a')
		print modelName+'_path_'+str(stepNum)+'.pdb'
		sys.stdout = fout

		os.chdir(currentPath)
	
		for sample in pathSeq:
			print "MODEL "+str(i)
			f = open(outputDir+"/"+modelName+"_new_"+str(sample)+".pdb", "r")
			print f.read()
			print "ENDMDL"
			i+=stepNum

		sys.stdout = orig_stdout
		fout.close()

		# Change to comboAnalsis
		currentPath=os.getcwd()
		d="comboAnalysis"
		if not os.path.exists(d):
		    os.makedirs(d)
		os.chdir(d)

		fout = file(modelName+'_combinedPath_'+str(stepNum)+'.pdb', 'a')
		sys.stdout = fout

		os.chdir(currentPath)
	
		for sample in pathSeq:
			print "MODEL "+str(i)
			f = open(outputDir+"/"+modelName+"_new_"+str(sample)+".pdb", "r")
			print f.read()
			print "ENDMDL"
			i+=stepNum

		sys.stdout = orig_stdout
		fout.close()

		# Change to comboAnalsis
		currentPath=os.getcwd()
		d="comboAnalysis"
		if not os.path.exists(d):
		    os.makedirs(d)
		os.chdir(d)

		fout = file(modelName+'_reversePath_'+str(stepNum)+'.pdb', 'a')
		sys.stdout = fout

		os.chdir(currentPath)
		
		for sample in reversePathSeq:
			print "MODEL "+str(ir)
			f = open(outputDir+"/"+targetName+"_new_"+str(sample)+".pdb", "r")
			print f.read()
			print "ENDMDL"
			ir+=stepNum

		sys.stdout = orig_stdout
		fout.close()

		#Finally, add the reverse path to the combined path
		reversePathList.reverse()
		reversePathSeq = reversePathList[::stepNum]


		# Change to comboAnalsis
		currentPath=os.getcwd()
		d="comboAnalysis"
		if not os.path.exists(d):
		    os.makedirs(d)
		os.chdir(d)

		fout = file(modelName+'_combinedPath_'+str(stepNum)+'.pdb', 'a')
		sys.stdout = fout

		os.chdir(currentPath)
	
		for sample in reversePathSeq:
			print "MODEL "+str(ir)
			f = open(outputDir+"/"+targetName+"_new_"+str(sample)+".pdb", "r")
			print f.read()
			print "ENDMDL"
			ir+=stepNum

		sys.stdout = orig_stdout
		fout.close()
	
	return 

def main():
	"""
	Checks arguments and writes the tree-path of a kgs pdb-file
	"""

	if len(sys.argv)<3:
		print "Usage: "+sys.argv[0]+" <stepNum> <path.pdb files in a row>"
		print "Use path files and only extract every <stepNum> configuration to write the path as a multi-model PDB."
		sys.exit(1)
		
	pathList=[]
	for i in range(len(sys.argv)-2):
		pathList.append(sys.argv[i+2])

	# writePath(sys.argv[1], sys.argv[2] if len(sys.argv)==3 else None)
	writePath(pathList, sys.argv[1])

if __name__ == "__main__":
	main()
