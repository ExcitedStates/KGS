#!/usr/bin/python
import pdb_structure
import sys

def getSequence(pdbFile):
	struc = pdb_structure.PDBFile(pdbFile) 
	print struc.getSequence()
	


def main():

	if len(sys.argv)!=2:
		print "Usage: "+sys.argv[0]+" <pdb-file>"
		sys.exit(1)
	
	getSequence(sys.argv[1])

if __name__ == "__main__":
        main()
