#!/usr/bin/python

import sys
import os
import numpy as np
import fnmatch
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from pylab import rcParams
rcParams['figure.figsize'] = 6.5, 6.5*3/4

ftSize=24
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

addLim = 0.3
forwardFinal=[]
revFinal=[]

if len(sys.argv)==1:
    print "Usage: pairwiseTreeDistances.py ", sys.argv[0]," <experiment directories in a row>"
    sys.exit(1)

numRuns = len(sys.argv)-1
run=0
#clashFac=75 #clashFac=[60, 75, 90]

allClashFreeDofsForward=[]
allClashFreeDofsReverse=[]

allDistancesForward=[]
allDistancesReverse=[]

currentPath = os.getcwd()
maxDofsFwd=0.0
minDofsFwd=99999
maxDofsRev=0.0
minDofsRev=99999

minDofs=99999
maxDofs = 0.0

fwdYLimit=9999
fwdXLimit=0.0
fwdMaxDist=0.0
revYLimit=0
revXLimit=99999
revMaxDist=0.0
iniDistance=1.0

relevant_path = sys.argv[1]
os.chdir(relevant_path)
regPlanner=1
dofs=0
countIt=0
minDistance = 99999

distancesForward=[]
distancesReverse=[]
clashFreeDofsForward=[]
clashFreeDofsReverse=[]

sample_ids = []
rmsdIni=[]
rmsdTarget=[]
reverse_sample_ids = []
reverse_rmsdIni=[]
reverse_rmsdTarget=[]

reverseIndex=0
forwardIndex=0

forwardName="init"
reverseName="target"

with open("output.txt") as outputFile:
    for line in outputFile:
        if "--initial " in line:
            forwardName = line[line.rfind("/")+1:line.rfind(".pdb")]
        if "--target " in line:
            reverseName = line[line.rfind("/")+1:line.rfind(".pdb")]
            
        if "Total DOFs:" in line:
            tokens = line.split(' ')
            totalDof= tokens[2]
            cycleDof = tokens[5]
            maxDofsFwd = int(totalDof[:-1]) - int(cycleDof[:-1])+int(tokens[9])
        if "Total DOFs in target:" in line:
            tokens = line.split(' ')
            totalDof= tokens[4]
            cycleDof = tokens[7]
            maxDofsRev = int(totalDof[:-1]) - int(cycleDof[:-1])+int(tokens[11])
        if "Initial Distance:" in line:
            tokens = line.split(' ')
            iniDistance = float(tokens[2])
            minDistance = iniDistance
            fwdMaxDist = iniDistance
            revMaxDist = iniDistance
 
        if "New structure:" in line:
            tokens = line.split(' ')
            # print tokens
            dofs =int(tokens[20])
            fileName=tokens[3]
            sampleID = int(fileName[fileName.rfind("_")+1:fileName.rfind(".pdb")])
            if( fileName.startswith(forwardName) ):
                clashFreeDofsForward.append(dofs)
                sample_ids.append(sampleID)
                distancesForward.append(minDistance/float(iniDistance))
                rmsd_ini=float(tokens[7])
                rmsdIni.append(rmsd_ini)
                if(rmsd_ini > fwdXLimit ):
                    fwdXLimit = rmsd_ini
                rmsd_target=float(tokens[11])
                rmsdTarget.append(rmsd_target)
                if(rmsd_target < fwdYLimit ):
                    fwdYLimit = rmsd_target
            else:
                clashFreeDofsReverse.append(dofs)
                reverse_sample_ids.append(sampleID)
                distancesReverse.append(minDistance/float(iniDistance))
                rmsd_ini=float(tokens[7])
                reverse_rmsdIni.append(rmsd_ini)
                if(rmsd_ini < revXLimit ):
                    revXLimit = rmsd_ini
                rmsd_target=float(tokens[11])
                reverse_rmsdTarget.append(rmsd_target)
                if(rmsd_target > revYLimit ):
                    revYLimit = rmsd_target
            if( dofs < minDofs):
                minDofs=dofs
            if( dofs > maxDofs):
                maxDofs=dofs
                
        if "Current shortest distance" in line:
            tokens = line.split(' ')
            distance = float(tokens[12])
            if distance < minDistance:
                minDistance=distance
            bestConf1 = int(tokens[8])
            bestConf2 = tokens[10]
            bestConf2 = int(bestConf2[:-1])
        
        if "The forward trajectory ends at" in line:
            tokens = line.split(' ')
            num = int(tokens[6])
            forwardFinal = num
            forwardIndex = sample_ids.index(forwardFinal)
            if num ==bestConf1:
                revFinal = bestConf2
            else:
                revFinal = bestConf1
            reverseIndex = reverse_sample_ids.index(revFinal)
                
allClashFreeDofsForward.append(list(clashFreeDofsForward))
allClashFreeDofsReverse.append(list(clashFreeDofsReverse))

allDistancesForward.append(list(distancesForward))
allDistancesReverse.append(list(distancesReverse))
os.chdir(currentPath)

print "Max dofs: ",maxDofs, ", min dofs: ",minDofs

currentPath = os.getcwd()
os.chdir(relevant_path)

included_extensions = ['*[0-9].pdb'] ;
outputPath = "output";

fig = plt.figure()
ax = fig.add_axes([0.11,0.17,0.75,0.75])

ax.set_xlim([0,revMaxDist+addLim])
ax.set_ylim([0,fwdMaxDist+addLim])

majorXLocator = MultipleLocator(1)
majorXFormatter = FormatStrFormatter('%d')
# minorXLocator = MultipleLocator(1)

majorYLocator = MultipleLocator(1)
majorYFormatter = FormatStrFormatter('%d')
# minorYLocator = MultipleLocator(1)

ax.xaxis.set_major_locator(majorXLocator)
ax.xaxis.set_major_formatter(majorXFormatter)
# ax.xaxis.set_minor_locator(minorXLocator)

ax.yaxis.set_major_locator(majorYLocator)
ax.yaxis.set_major_formatter(majorYFormatter)

ax.set_ylabel("RMSD to target [\AA]",fontsize=ftSize)
ax.set_xlabel("RMSD to initial [\AA]",fontsize=ftSize)

jet = cm = plt.get_cmap('jet') 
#cNorm  = colors.Normalize(vmin=minDofs, vmax=maxDofs)
cNorm  = colors.Normalize(vmin=0, vmax=maxDofs-minDofs)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#print scalarMap.get_clim()

dots=[]
clashFreeDofsForward=allClashFreeDofsForward[run]
clashFreeDofsReverse=allClashFreeDofsReverse[run]
distancesForward=allDistancesForward[run]
distancesReverse=allDistancesReverse[run]

for i in range(len(clashFreeDofsForward)):
    colorval=scalarMap.to_rgba(maxDofs-clashFreeDofsForward[i])
    #colorval=scalarMap.to_rgba(clashFreeDofsForward[i])
    dot, = ax.plot(rmsdIni[i],rmsdTarget[i],'v',color=colorval,markeredgecolor=colorval,alpha=0.7)
    dots.append(dot)
    
for i in range(len(clashFreeDofsReverse)):
    colorval=scalarMap.to_rgba(maxDofs-clashFreeDofsReverse[i])
    #colorval=scalarMap.to_rgba(clashFreeDofsReverse[i])
    dot, = ax.plot(reverse_rmsdTarget[i],reverse_rmsdIni[i],'^',color=colorval,markeredgecolor=colorval,alpha=0.7)
    dots.append(dot)
 
colorval=scalarMap.to_rgba(maxDofs-clashFreeDofsForward[forwardIndex])
dot, = ax.plot(rmsdIni[forwardIndex],rmsdTarget[forwardIndex],'v',color=colorval,zorder=3,)
colorval=scalarMap.to_rgba(maxDofs-clashFreeDofsReverse[reverseIndex])
dot, = ax.plot(reverse_rmsdTarget[reverseIndex],reverse_rmsdIni[reverseIndex],'^',color=colorval,zorder=3)

print "Final forward: "+str(forwardFinal)+", best reverse: "+str(revFinal)
    
scalarMap.set_array(dots)

os.chdir(currentPath)
d="comboAnalysis"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir(d)   

#fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.22, 0.02, 0.6])
cbar=fig.colorbar(scalarMap, cax=cbar_ax)
cbar.ax.tick_params(labelsize=ftSize)
cbar.set_ticks([0,maxDofs-minDofs])
cbar.set_ticklabels([0,maxDofs-minDofs])
# cbar.yaxis.set_major_formatter(majorYFormatter)
#cbar_ax.xaxis.set_ticks_position('top')
cbar.set_label('clash constraints',labelpad=-10,fontsize=ftSize)

# plt.savefig("pairedRMSD_twoTrees_"+str(clashFac)+".png", transparent='True',format='png',dpi=600,figsize=(4,4))
plt.savefig("pairedRMSD_twoTrees.png", transparent='True',format='png',dpi=600,figsize=(6.5, 6.5*3/4)) #figsize=(4,4)

os.chdir(currentPath)
