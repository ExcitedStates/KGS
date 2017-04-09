#!/usr/bin/python
"""
KGS generates PDB-files with a 'REMARK Tree-path' in the header that can be used to 
reconstruct the exploration path from any file to the initial. This program reads a PDB-file
parses the path and prints the contents of all PDB-files in the path (assumed to be 
"""

import sys

def writePath(forwardPDBFile, reversePDBFile, initialPath=None):
	pathStr=""
	reversePathStr=""
	forwardPathPrint=""
	reversePathPrint=""
	
	samplePreFwd = forwardPDBFile[0:forwardPDBFile.rfind("_")+1]
	samplePreRev = reversePDBFile[0:reversePDBFile.rfind("_")+1]

	fp = open(forwardPDBFile)
	for line in fp:
		if "REMARK	Tree-path" in line:
			pathStr = line[line.find("=")+2:].strip()
			forwardPathPrint=line
			break;

	fp.close()
	pathList = pathStr.split(" ")
	del pathList[-1] #Remove last element
	pathList.reverse()
	
	fp = open(reversePDBFile)
	for line in fp:
		if "REMARK	Tree-path" in line:
			reversePathStr = line[line.find("=")+2:].strip()
			reversePathPrint="REMARK	Reverse path follows "+line[line.find("=")+2:]
			break;

	fp.close()
	reversePathList = reversePathStr.split(" ")
	del reversePathList[-1] #Remove starting element

	i = 1
	if initialPath!=None:
		print "MODEL "+str(i)
		f = open(initialPath, "r")
		print f.read()
		print "ENDMDL"
		i+=1

	print forwardPathPrint #Print forward path information
	
	for sample in pathList:
		print "MODEL "+str(i)
		f = open(samplePreFwd+sample+".pdb", "r")
		print f.read()
		print "ENDMDL"
		i+=1

	print reversePathPrint #Print reverse path information
	
	for sample in reversePathList:
		print "MODEL "+str(i)
		f = open(samplePreRev+sample+".pdb", "r")
		print f.read()
		print "ENDMDL"
		i+=1

def main():
	"""
	Checks arguments and writes the tree-path of a kgs pdb-file to stdout
	"""

	if len(sys.argv)<2:
		print "Usage: "+sys.argv[0]+" <kgs forward pdb-file> <kgs reverse pdb-file> [<initial pdb-file>]"
		print "If the 'REMARK Tree-Path' line is found a pdb-file, the entire path is written to screen as a multi-model PDB."
		sys.exit(1)

	writePath(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv)==4 else None)

if __name__ == "__main__":
	main()
