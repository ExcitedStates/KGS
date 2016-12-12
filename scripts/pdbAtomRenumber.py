#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  

__description__ = \
"""
pdbAtomRenumber.py

Renumbers atom records in a pdb without touching residue numbers.
"""

__author__ = "Dominik Budday"
__date__ = "161021"

import os, sys


"""
Renumber all atoms in pdb file, starting from 1.
"""


if len(sys.argv)< 2:
    print "Usage: pdbAtomRenumber.py ", sys.argv[0]," <inputPDBfiles in a row>"
    
# stfRes = []
# kgsRes = {}
resList = {}
noShowVal = 0
#List of residues: Strong = 2; Medium = 1; Weak = 0; not present
# conectMap={}
firstHet=0
firstConect=0
for input in sys.argv[1:]:
    counter=1
    with open(input) as pdb:
        output = input[:-4]+"_renum.pdb"
        outputFile = open(output, 'w')
        for line in pdb:
            # For ATOM and HETATM record, update residue number
            if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
                line = "%s%5s%s" % (line[0:6],counter,line[11:])
                counter += 1
            if line[0:6] == "HETATM":
                # oldID = int(line[6:11])
                # conectMap[oldID]=counter
                if firstHet==0:
                    firstHet = counter
                line = "%s%5s%s" % (line[0:6],counter,line[11:])
                counter += 1
            #Update CONECT records
            if line[0:6] == "CONECT":
                newLine = "CONECT "
                tokens=line.split(' ')
                for entry in tokens[1:]:
                    if firstConect == 0:
                        firstConect = int(entry)
                        diff = firstConect - firstHet
                    # newID = conectMap[int(entry)]
                    newID = int(entry) - diff
                    newLine = newLine+str(newID)+" "
                line = newLine+"\n"
            outputFile.write(line)
        outputFile.close()