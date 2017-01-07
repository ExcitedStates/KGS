#!/usr/bin/python
"""
KGS generates PDB-files with a 'REMARK Tree-path' in the header that can be used to 
reconstruct the exploration path from any file to the initial. This program reads a PDB-file
parses the path and prints the contents of all PDB-files in the path (assumed to be 
"""

import sys

def writePath(pdbPath, initialPath=None):
	pathFileSepIdx = pdbPath.rfind("/")
	#outputDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
	sampleNum = int(pdbPath[pdbPath.rfind("_")+1:pdbPath.rfind(".pdb")])

	samplePre = pdbPath[0:pdbPath.rfind("_")+1]

	fp = open(pdbPath)
	for line in fp:
		if "REMARK	Tree-path" in line:
			pathStr = line[line.find("=")+2:].strip()
			break;

	fp.close()
	pathList = pathStr.split(" ")
	del pathList[-1] #Remove last element
	#pathList = map(int, pathList) #Convert to list-of-int
	#pathList = sorted(pathList) #Reverse
	pathList.reverse()

	i = 1
	if initialPath!=None:
		print "MODEL "+str(i)
		f = open(initialPath, "r")
		print f.read()
		print "ENDMDL"
		i+=1

	for sample in pathList:
		print "MODEL "+str(i)

		#f = open(outputDir+"/newpdb_"+sample+".pdb", "r")
		f = open(samplePre+sample+".pdb", "r")
		print f.read()
		print "ENDMDL"
		i+=1


def main():
	"""
	Checks arguments and writes the tree-path of a kgs pdb-file to stdout
	"""

	if len(sys.argv)<2:
		print "Usage: "+sys.argv[0]+" <kgs pdb-file> [<initial pdb-file>]"
		print "If the 'REMARK Tree-Path' line is found a pdb-file, the entire path is written to screen as a multi-model PDB."
		sys.exit(1)

	writePath(sys.argv[1], sys.argv[2] if len(sys.argv)==3 else None)

if __name__ == "__main__":
	main()
