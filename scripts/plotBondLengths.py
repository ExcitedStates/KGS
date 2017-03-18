#!/usr/bin/python

import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

if len(sys.argv)< 5:
    print "Usage: ", sys.argv[0]," <_ini_allCovBonds.txt file> <_target_allCovBonds.txt file> <_ini_allAnchors.txt file> <_target_allAnchors.txt file> <constraintLog> (optional)"

def plotLengths(bonds,bondLengthsDiff,saveFileName=None):
    import matplotlib.pyplot as plt
    
    plt.plot(bonds,bondLengthsDiff,'b-')
    plt.xlabel('# of bonds')
    plt.ylabel('bond length difference')
    if saveFileName:
        plt.savefig(saveFileName)
    else:
        plt.show()
    plt.clf()

bondLengthsIni = [];
bondLengthsTarget = [];
bondLengthsDiff=[]

path_file_ini = sys.argv[1];
path_file_target = sys.argv[2];

print path_file_ini
print path_file_target

#Get bond lengths ini
with open(path_file_ini) as iniFile:
    for line in iniFile:
        bondLengthsIni.append(float(line))

#Get bond lengths target
with open(path_file_target) as targetFile:
    for line in targetFile:
        bondLengthsTarget.append(float(line))

if( len(bondLengthsIni) != len(bondLengthsTarget)):
    print "Not the same number of bonds, quitting!"
    sys.exit()
    
for index in range(len(bondLengthsIni)):
    bondLengthsDiff.append(abs(bondLengthsIni[index]-bondLengthsTarget[index]))

bonds = np.arange(len(bondLengthsIni))

currentPath = os.getcwd()
d="constraintPlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("constraintPlots")

plotLengths(bonds,bondLengthsDiff,saveFileName="covBondLengths.png")

os.chdir(currentPath)


## Do the same for the hydrogen bonds

bondLengthsIni = [];
bondLengthsTarget = [];
bondLengthsDiff=[]

path_file_ini = sys.argv[3];
path_file_target = sys.argv[4];

#Get bond lengths ini
with open(path_file_ini) as iniFile:
    for line in iniFile:
        bondLengthsIni.append(float(line))

#Get bond lengths target
with open(path_file_target) as targetFile:
    for line in targetFile:
        bondLengthsTarget.append(float(line))

if( len(bondLengthsIni) != len(bondLengthsTarget)):
    print "Not the same number of hydrogen bonds, quitting!"
    sys.exit()

for index in range(len(bondLengthsIni)):
    bondLengthsDiff.append(abs(bondLengthsIni[index]-bondLengthsTarget[index]))

bonds = np.arange(len(bondLengthsIni))

bondIds=[]
if(len(sys.argv)==6):
    constraint_file_name = sys.argv[5]
    print constraint_file_name

    with open(constraint_file_name) as constraintFile:
        for line in constraintFile:
            if line == '\n':
                break;
            tokens = line.split(' ')
            bondId = int(tokens[3])
            id1= int(tokens[5])
            id2 = int(tokens[7])
            bondIds.append([id1,id2])


currentPath = os.getcwd()
d="constraintPlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("constraintPlots")

plotLengths(bonds,bondLengthsDiff,saveFileName="anchorLengths.png")

for i in range(len(bondLengthsDiff)):
    if bondLengthsDiff[i] > 0.2:
        print str(bondLengthsDiff[i])+" AA length difference at bond "+str(bonds[i])+" between atoms "+str(bondIds[i])

