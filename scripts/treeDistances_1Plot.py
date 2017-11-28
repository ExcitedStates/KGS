#!/usr/bin/python

import sys
import os
import numpy as np
import fnmatch
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.pyplot import *

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from pylab import rcParams
rcParams['figure.figsize'] = 6.5, 4.5

ftSize=12
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

# def plotDistance(fig,samples,distances, legendName,saveFileName=None):
#     import matplotlib.pyplot as plt
#     
#     # plt.plot(samples,distances, 'b-')
#     plt.plot(samples,distances, label=legendName)
#     plt.plot(samples[0::1000],distances[0::1000], 'rs')
#     plt.xlabel("# samples in 1000")
#     plt.ylabel("RMSD between trees")
#     plt.xticks(np.arange(0,len(samples)+1,1000),['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
#     ax = gca;
#     ax.XTickLabel = {'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'};
#     plt.grid(which='major',axis='both')
#     plt.legend()
#     
#     if saveFileName:
#         plt.savefig(saveFileName, transparent='True',format='png',dpi=600)
#     else:
#         plt.show()
#     plt.show()
#     # plt.clf()

if len(sys.argv)<6:
    print "Usage: treeDistances.py ", sys.argv[0]," <number Proteins> <number rounds>, < <protein name> <output directory> <numSamples in 1000> > repeat"
    sys.exit()
    
maxLength=0
clashFac=[0.60,0.75,0.90]
#colorList=['b','r','g','y','c','m','k']
colorList=['red','green','blue','yellow','cyan','magenta','black','orange','darkolivegreen','chocolate']
symbolList=['s','o','d','v','^','<','>','*','.',',']

numProteins = int(sys.argv[1])
numRounds = int(sys.argv[2])
remainingInput = sys.argv[3:]
print numProteins," Proteins, ",numRounds," rounds"
plots=[0]*numProteins

#fig, axarr = plt.subplots((numRounds+1)/2, 2)
#fig, axarr = plt.subplots(1,numRounds)
fig = plt.figure()

for rounds in range(numRounds):
    #ax = axarr[round(rounds/2),rounds%2]
    #ax = axarr[rounds]
    #ax = subplot((numRounds+1)/2,2,rounds+1)
    if rounds==0:
        ax = fig.add_axes([0.09,0.59,0.4,0.38])
    if rounds==1:
        ax = fig.add_axes([0.58,0.59,0.4,0.38])
    if rounds==2:
        ax = fig.add_axes([0.09,0.1,0.4,0.38])
    
    for inp in range(numProteins):
        proteinName=remainingInput[inp*3]
        relevant_path = remainingInput[inp*3+1]
        baseName = relevant_path[:-1]
        numSamples = int(remainingInput[inp*3+2])*1000 #given index number in thousand samples, e.g. 10 --> 10000
    
        distances=[]
        samples=np.arange(0,numSamples+1)
        
        currentPath = os.getcwd()
        #print currentPath
        os.chdir(relevant_path)
        
        count=0
        iniDistance=1.0
     
        os.chdir(currentPath)
        # os.chdir("..")
        iniDistance=1.0
        minDistance = 99999
        newPath = baseName
        #print newPath
        os.chdir(newPath)
        count=0
        with open("output.txt") as outputFile:
            for line in outputFile:
                if "Current shortest distance" in line:
                    tokens = line.split(' ')
                    distanceVal = tokens[12]
                    distance =float(distanceVal[:-1])
                    if count ==0:
                        iniDistance=distance
                        count=1
                    if distance < minDistance:
                        distances.append(distance/float(iniDistance))
                        minDistance=distance
                    else:
                        distances.append(minDistance/float(iniDistance))
            print iniDistance, minDistance, minDistance/float(iniDistance)
                        
        os.chdir(currentPath)
        # os.chdir("..") 
        
        # print len(samples), len(distances)
        while( len(samples) > len(distances)):
            distances.append(distances[-1])
        maxLength=max(maxLength,len(samples))
        os.chdir(currentPath)
        # Option with markers enabled
        # plots[inp], =ax.plot(samples,distances, markevery=1000,marker=symbolList[inp],markeredgecolor=colorList[inp],color=colorList[inp],label=proteinName, linewidth=2)
        plots[inp], =ax.plot(samples,distances, color=colorList[inp],label=proteinName, linewidth=2)
        
    majorTicks = np.arange(0,maxLength,5000)
    majorTickLabels = np.arange(0,(maxLength)/1000+1,5)
    #print maxLength/5000
    #print len(samples),len(distances),maxLength
    #print majorTicks
    #print majorTickLabels
    majorYTicks = np.arange(0,1.1,0.2)
    ax.set_xticks(majorTicks)
    ax.set_xticklabels(majorTickLabels)
    ax.set_yticks(majorYTicks)
    ax.set_xlabel("\# samples in 1000",fontsize=ftSize)
    ax.set_ylabel("fRMSD between trees",fontsize=ftSize)
    ax.grid(axis='both')
    ax.text(15100,0.85,'$c_f=$'+str(clashFac[rounds]),fontsize=ftSize)
    remainingInput=remainingInput[numProteins*3:]
        
#axarr[(numRounds+1)/2-1, 1].legend(loc=7,fontsize=ftSize)
# axarr[numRounds].legend(loc=7,fontsize=ftSize)

legend(bbox_to_anchor=(1.25, 0.05, 0.95, 0.9), loc=5,
       ncol=2, mode="expand", borderaxespad=0.)

remainingInput = sys.argv[3:3+numProteins*3]
#l1=axarr[0].legend([plots[0],(plots[0],plots[1]),(plots[0],plots[1],plots[2]),(plots[0],plots[1],plots[2],plots[3])],[remainingInput[0],remainingInput[3],remainingInput[6],remainingInput[9]],loc=7,fontsize=ftSize)
##l2=axarr[1].legend([plots[3],(plots[3],plots[4]),(plots[3],plots[4],plots[5])],[remainingInput[12],remainingInput[15],remainingInput[18]],loc=7,fontsize=ftSize)


d="comboAnalysis"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("comboAnalysis")   
saveFileName="rmsdBetweenTrees"
saveFormat='png'

if saveFileName:
    # plt.savefig(saveFileName+'.'+saveFormat, transparent='True',format=saveFormat,dpi=600)
    plt.savefig(saveFileName+'.'+saveFormat,format=saveFormat,dpi=600)
else:
    plt.show()
plt.show()
    

                