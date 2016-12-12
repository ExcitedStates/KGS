#!/usr/bin/python
import sys

pdbs = sys.argv[1:len(sys.argv)]
for i in range(0,len(pdbs)):
	print "MODEL "+str(i)
	f = open(pdbs[i], "r")
	print f.read()
	print "ENDMDL"

