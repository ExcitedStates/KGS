#!/usr/bin/python

import sys
import os
import numpy as np
import fnmatch
import operator
import matplotlib as mpl
mpl.use('Agg')
from combineKGSPath_steps import extractPath

if len(sys.argv)==1:
    print "Usage: ", sys.argv[0]," <path.pdb file>"

meanStrains = []
maxStrains = []
reverse_meanStrains = []
reverse_maxStrains = []

bins=[0, 0.1, 1, 10, 100]

path_file = sys.argv[1];
#relevant_path = sys.argv[2]

constraint_file_name = "kgs_report.log"
outputFile = "output.txt"

#Get max distance constraint violation and path length
# with open(path_file) as pdbPath_file:
#     for line in pdbPath_file:
#         if "REMARK\tFound connected" in line:
#             tokens = line.split(' ')
#             pathLength = int(tokens[4])+1
#             print "Path length: ",pathLength
        
#         if "REMARK\tPath is" in line:
#             tokens = line.split(' ')
#             for entry in reversed(tokens[2:]):
#                 path_ids.append(int(entry[:-1]))

path_ids, reverse_path_ids = extractPath(path_file)

pathLength = len(path_ids)
reverse_pathLength = len(reverse_path_ids)

#Get the number of hbonds
fwdHbonds=True
numHBonds=0
targetNumHBonds=0
with open(outputFile) as log_file:
    for line in log_file:
        if "Molecule has:" in line:
            fwdHbonds=True
            continue
        if "Target has:" in line:
            fwdHbonds=False
            continue
        if "hydrogen bonds" in line:
            tokens = line.split(' ')
            if fwdHbonds == True:
                numHBonds = int(tokens[1])
                print "Number of hbonds: ",numHBonds
            else:
                targetNumHBonds = int(tokens[1])
                print "Number of target hbonds: ",targetNumHBonds               
                break

if(pathLength < 2):
    print "No path encountered, returning."
    sys.exit() 

meanStrains.append(0)
maxStrains.append(0)
pathLengths=np.arange(pathLength)

currentConf = -1;
currentPathEntry = 1;
meanStrain=0
maxStrain=0
binData=[0]*numHBonds
strainedBonds=[0, 0]*numHBonds
secStrainedBonds=[0, 0]*numHBonds
allBonds = {}
foundStrain=0

with open(constraint_file_name) as constraintFile:
    for line in constraintFile:
        if line.startswith("Conformation"):
            currentConf = int(line[13:])
            continue
        if(currentConf != path_ids[currentPathEntry-1]):
            continue
        tokens = line.split(' ')
        bondId = int(tokens[3])
        id1= int(tokens[5])
        id2 = int(tokens[7])
        
#        print "At configuration: ",currentConf
#        print "At hbond: ",bondId
        strain = float(tokens[11])
        meanStrain += abs(strain)
        if(abs(strain) > maxStrain):
            maxStrain = abs(strain)
        if(abs(strain)>binData[bondId-1]):
            binData[bondId-1]=abs(strain)
        # if(abs(strain)>bins[-4]):
        #     print "Bin 0.1 to 1: ",id1,",",id2
        if(abs(strain)>bins[-3]):
            secStrainedBonds[bondId]=[id1, id2]
            foundSec=1
        if(abs(strain)>bins[-2]):
            strainedBonds[bondId]=[id1, id2]
            foundStrain=1

        # all bonds in order
        cons = (id1,id2)
        if cons in allBonds:
        	oldVal = allBonds[cons]
        	if ( abs(strain) > oldVal ):
        		allBonds[cons] = abs(strain)
        else:
        	allBonds[cons] = abs(strain)
            
        if(bondId == numHBonds):#reached end of current conf
            meanStrains.append(meanStrain/numHBonds)
            maxStrains.append(maxStrain)
            maxStrain=0
            meanStrain=0
            currentPathEntry += 1
        if currentPathEntry == pathLength:
            break
print "Forward path"
# if foundStrain==1:
#     print "Bin > 10 %"
#     print strainedBonds
# else:
#     print "Bin 1 to 10 %"
#     print secStrainedBonds

sorted_bondStrains = sorted(allBonds.items(),key=operator.itemgetter(1))
sorted_bondStrains.reverse()
print "Top ten strained bonds"
for i in range(0,min(len(sorted_bondStrains),10)):
    entry =  sorted_bondStrains[i]
    key = entry[0]
    print "distance strained = id "+str(key[0])+", id "+str(key[1])
    if key[0] ==1382:
        print "Strain at id "+str(key[0])+", id "+str(key[1])+": "+str(entry)


#Same for revese sampling
reverse_meanStrains.append(0)
reverse_maxStrains.append(0)
reverse_pathLengths=np.arange(reverse_pathLength)

currentConf = -1;
currentPathEntry = 1;
reverse_meanStrain=0
reverse_maxStrain=0
reverse_binData=[0]*targetNumHBonds
reverse_strainedBonds=[0, 0]*targetNumHBonds
reverse_secStrainedBonds=[0, 0]*targetNumHBonds
reverse_allBonds = {}
foundStrain=0
print reverse_path_ids
with open(constraint_file_name) as constraintFile:
    for line in constraintFile:
        if line.startswith("Conformation"):
            currentConf = int(line[13:])
            continue
        if(currentConf != reverse_path_ids [currentPathEntry-1]):
            continue
        tokens = line.split(' ')
        bondId = int(tokens[3])
        id1= int(tokens[5])
        id2 = int(tokens[7])
        
#        print "At configuration: ",currentConf
#        print "At hbond: ",bondId
        strain = float(tokens[11])
        reverse_meanStrain += abs(strain)
        if(abs(strain) > reverse_maxStrain):
            reverse_maxStrain = abs(strain)
        if(abs(strain)>reverse_binData[bondId-1]):
            reverse_binData[bondId-1]=abs(strain)
        # if(abs(strain)>bins[-4]):
        #     print "Bin 0.1 to 1: ",id1,",",id2
        if(abs(strain)>bins[-3]):
            reverse_secStrainedBonds[bondId]=[id1, id2]
            foundSec=1
        if(abs(strain)>bins[-2]):
            reverse_strainedBonds[bondId]=[id1, id2]
            foundStrain=1

        # all bonds in order
        cons = (id1,id2)
        if cons in reverse_allBonds:
            oldVal = reverse_allBonds[cons]
            if ( abs(strain) > oldVal ):
                reverse_allBonds[cons] = abs(strain)
        else:
            reverse_allBonds[cons] = abs(strain)
            
        if(bondId == targetNumHBonds):#reached end of current conf
            reverse_meanStrains.append(reverse_meanStrain/targetNumHBonds)
            reverse_maxStrains.append(reverse_maxStrain)
            reverse_maxStrain=0
            reverse_meanStrain=0
            currentPathEntry += 1
        if currentPathEntry == reverse_pathLength:
            break

print "REVERSE PATH"
# if foundStrain==1:
#     print "Bin > 10 %"
#     print reverse_strainedBonds
# else:
#     print "Bin 1 to 10 %"
#     print reverse_secStrainedBonds

sorted_reverse_bondStrains = sorted(reverse_allBonds.items(),key=operator.itemgetter(1))
sorted_reverse_bondStrains.reverse()
print "Top ten strained bonds"
for i in range(0,min(len(sorted_reverse_bondStrains),10)):
    entry =  sorted_reverse_bondStrains[i]
    key = entry[0]
    print "distance strained = id "+str(key[0])+", id "+str(key[1])
    if key[0] ==1382:
        print "Strain at id "+str(key[0])+", id "+str(key[1])+": "+str(entry)
    
def plotConstraintViolation(xData,meanY,maxY,saveFileName=None):
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    ax1.plot(xData,meanY,'b-')
    ax1.set_xlabel('path length in # of samples')
    ax1.set_ylabel('mean absolute distance violation in %', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax2 = ax1.twinx()
    ax2.plot(xData,maxY,'r-')
    ax2.set_ylabel('max absolute distance violation in %', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    if saveFileName:
        plt.savefig(saveFileName)
    else:
        plt.show()
        
def plotHistogram(binData,bins,saveFileName=None):
    import matplotlib.pyplot as plt
    import math
    histArray=np.histogram(binData,bins)
    yData=histArray[0]
    #yData.append(yData[-1])
    numBins=len(bins)
    xVals=np.arange(len(bins)-1)
    print xVals
    print histArray[0]
    fig, ax1 = plt.subplots()
    ax1.step(xVals,yData,where='post')
    i=0
    majorTicks = np.arange(numBins)
    majorTickLabels = bins
    print majorTickLabels
    #print maxLength/1000
    #print len(samples),len(distances),maxLength
    #print majorTicks
    #print majorTickLabels
    ax1.set_xticks(majorTicks)
    ax1.set_xticklabels(majorTickLabels)
    ax1.set_xlabel('maximum distance violation in %')
    ax1.set_ylabel('number of bonds')
    if saveFileName:
        plt.savefig(saveFileName)
    else:
        plt.show()

currentPath = os.getcwd()
d="constraintPlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("constraintPlots")

#print" Path: ",pathLengths
#print" Means: ",meanStrains
#print "Maxs: ",maxStrains

plotConstraintViolation(pathLengths,meanStrains,maxStrains,saveFileName="forwardPath_meanMaxViolation.png")

#plotHistogram(binData,bins,saveFileName="hist_maxBondViolation.png")


plotConstraintViolation(reverse_pathLengths,reverse_meanStrains,reverse_maxStrains,saveFileName="reversePath_meanMaxViolation.png")

#plotHistogram(reverse_binData,bins,saveFileName="histRev_maxBondViolation.png")


os.chdir(currentPath)