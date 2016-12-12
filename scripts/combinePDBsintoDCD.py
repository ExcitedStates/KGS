#!/usr/bin/env python
# USAGE: combinePDBsintoDCD.py out.dcd `ls -rth *pdb`
import sys
from subprocess import call

catdcd = "catdcd"   # path to catdcd exec 
trajout = sys.argv[1]
pdbs = sys.argv[2:len(sys.argv)]
s = ""
for i in range(0,len(pdbs)):
	s += " "+pdbs[i]

s = catdcd + " -o " + trajout + " -pdb " + s + " > /dev/null"
call( s, shell=True )

