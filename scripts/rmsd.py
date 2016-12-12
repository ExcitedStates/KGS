#!/usr/bin/python

import sys
from pdb_structure import *

if __name__ == '__main__':
	names = ["C4'","CA"]
	userNames = []
	for a in range(len(sys.argv)-1,0,-1):
		if sys.argv[a]=="-name":
			userNames.append(sys.argv[a+1])
			del sys.argv[a]
			del sys.argv[a]
	if userNames: names = userNames

	if len(sys.argv)==3:
		s1 = PDBFile(sys.argv[1])
		s2 = PDBFile(sys.argv[2])
		print s1.rmsd(s2, names=names)
	elif len(sys.argv)>=3:
		s1 = PDBFile(sys.argv[1])
		for i in range(2,len(sys.argv)):
			s2 = PDBFile(sys.argv[i])
			print sys.argv[i],":",s1.rmsd(s2, names=names)
	else:
		print "Usage:",sys.argv[0],"[-name <atomname>] <pdb-file1> <pdb-file2> [...]"

