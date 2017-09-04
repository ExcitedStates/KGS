#!/usr/bin/python
import sys
import os.path

pdbs = sys.argv[1:len(sys.argv)]
i = 0
for pdbfname in pdbs:
    try:
        f = open(pdbfname, "r")
        print "MODEL "+str(i)
        print f.read()
        print "ENDMDL"
        i += 1
    except IOError:
        pass

