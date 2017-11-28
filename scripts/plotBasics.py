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
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

if len(sys.argv)==1:
    print "Usage: plotDistances.py ", sys.argv[0]," <output directory>"

sample_ids = []
distanceToIni = []
distanceToParent =[]
distanceToTarget = []
energy=[]
clashFreeDofs=[]

reverse_sample_ids = []
reverse_distanceToIni = []
reverse_distanceToParent =[]
reverse_distanceToTarget = []
reverse_energy=[]
reverse_clashFreeDofs=[]

relevant_path = sys.argv[1]
included_extensions = ['*[0-9].pdb'] ;


regPlanner=1
maxDofs=0
minDofs=99999
dofs=0
forwardName="init"
targetName="target"
with open("output.txt") as outputFile:
    for line in outputFile:
        if "--initial " in line:
            forwardName = line[line.rfind("/")+1:line.rfind(".pdb")]
        if "--target " in line:
            targetName = line[line.rfind("/")+1:line.rfind(".pdb")]
        if "Initial Distance:" in line:
            tokens = line.split(' ')
            initialDistance = float(tokens[2])
            #Add initial
            # sample_ids.append(0)
            # distanceToIni.append(0.0)
            # distanceToParent.append(0.0)
            # distanceToTarget.append(iniDistance)
            break;

all_files =[fn for fn in os.listdir(relevant_path) if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) ];
file_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and fn.startswith(forwardName)];
target_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and fn.startswith(targetName)];

currentPath = os.getcwd()
os.chdir(relevant_path)

for pdb_file_name in file_names:
        sample_id = -1
        sample_parent = -1
        iniDistance = 0.0
        parentDistance = 0.0
        targetDistance = 1000.0
        clashFreeDof = 0
        en = 0
        with open(pdb_file_name) as pdb_file:
            for line in pdb_file:
                if "REMARK\tID" in line:
                    sample_id = int(line[12:])
                if "REMARK\tParent ID" in line:
                    sample_parent = int(line[19:])
                if "REMARK\tDistance_initial" in line:
                    iniDistance = float(line[26:])
                if "REMARK\tDistance to parent" in line:
                    parentDistance = float(line[28:])
                if "REMARK\tDistance to target" in line:
                    targetDistance = float(line[28:])
                if "REMARK\tClash free dofs" in line:
                    clashFreeDof = float(line[25:])
                if "REMARK\tVdw energy" in line:
                    en = float(line[20:])
                    break

        sample_ids.append(sample_id)
        distanceToIni.append(iniDistance)
        distanceToTarget.append(targetDistance)
        distanceToParent.append(parentDistance)
        clashFreeDofs.append(clashFreeDof)
        energy.append(en)
        if clashFreeDof > maxDofs:
            maxDofs = clashFreeDof
        if clashFreeDof < minDofs:
            minDofs = clashFreeDof
        
for pdb_file_name in target_names:
#        print pdb_file_name
        reverse_sample_id = -1
        reverse_sample_parent = -1
        reverse_iniDistance = 0.0
        reverse_parentDistance = 0.0
        reverse_targetDistance = 1000.0
        reverse_clashFreeDof = 0
        reverse_en = 0
        with open(pdb_file_name) as pdb_file:
            for line in pdb_file:
                if "REMARK\tID" in line:
                    reverse_sample_id = int(line[12:])
                if "REMARK\tParent ID" in line:
                    reverse_sample_parent = int(line[19:])
                if "REMARK\tDistance_initial" in line:
                    reverse_iniDistance = float(line[26:])
                if "REMARK\tDistance to parent" in line:
                    reverse_parentDistance = float(line[28:])
                if "REMARK\tDistance to target" in line:
                    reverse_targetDistance = float(line[28:])
                if "REMARK\tClash free dofs" in line:
                    reverse_clashFreeDof = float(line[25:])
                if "REMARK\tVdw energy" in line:
                    reverse_en = float(line[20:])
                    break
        
        reverse_sample_ids.append(reverse_sample_id)
        reverse_distanceToIni.append(reverse_iniDistance)
        reverse_distanceToTarget.append(reverse_targetDistance)
        reverse_distanceToParent.append(reverse_parentDistance)
        reverse_clashFreeDofs.append(reverse_clashFreeDof)
        reverse_energy.append(reverse_en)
        if reverse_clashFreeDof > maxDofs:
            maxDofs = reverse_clashFreeDof
        if reverse_clashFreeDof < minDofs:
            minDofs = reverse_clashFreeDof
            
print "Max dofs: ",maxDofs, ", min dofs: ",minDofs

def plotDistance(xData,samples, yString,saveFileName=None):
    import matplotlib.pyplot as plt
    plt.plot(xData,samples, 'ro')
    plt.xlabel("# samples")
    plt.ylabel(yString)
    if saveFileName:
        plt.savefig(saveFileName, transparent='True')
    else:
        plt.show()
    plt.clf()


os.chdir("..")
d="distancePlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("distancePlots")

plotDistance(sample_ids,distanceToIni,"RMSD to initial",saveFileName="distanceToIni.png")
plotDistance(sample_ids,distanceToTarget,"RMSD to target",saveFileName="distanceToTarget.png")

plt.figure()
plotDistance(reverse_sample_ids,reverse_distanceToIni,"RMSD to reverse initial",saveFileName="reverse_distanceToIni.png")
plotDistance(reverse_sample_ids,reverse_distanceToTarget,"RMSD to reverse target",saveFileName="reverse_distanceToTarget.png")
  
fig = plt.figure()
ax = fig.add_axes([0.1,0.58,0.9,0.4])

jet = cm = plt.get_cmap('RdBu') 
cNorm  = colors.Normalize(vmin=minDofs, vmax=maxDofs)
# cNorm  = colors.Normalize(vmin=minDofs, vmax=maxDofs)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
#print scalarMap.get_clim()

dots=[]
for i in range(len(clashFreeDofs)):
    colorval=scalarMap.to_rgba(clashFreeDofs[i])
    dot, = ax.plot(distanceToIni[i],distanceToTarget[i],'s',color=colorval)
    dots.append(dot)
 
scalarMap.set_array(clashFreeDofs)
cbar = fig.colorbar(scalarMap)
#cbar = fig.colorbar(dots)
cbar.ax.set_yticks([-1,1])
cbar.update_ticks()
s=str(minDofs/float(maxDofs))
#cbar.set_label('mobility due to constraints')
cbar.set_label('remaining dofs')
#cbar.ax.set_yticklabels(['',''])
cbar.ax.set_yticklabels([str(minDofs),str(maxDofs)])
cbar.update_ticks()
    
plt.xlabel("RMSD to ini")
plt.ylabel("RMSD to target")

ax = fig.add_axes([0.1,0.09,0.9,0.4])
for i in range(len(reverse_distanceToTarget)):
    colorval=scalarMap.to_rgba(reverse_clashFreeDofs[i])
    dot, = ax.plot(reverse_distanceToTarget[i],reverse_distanceToIni[i],'s',color=colorval)
    dots.append(dot)

scalarMap.set_array(reverse_clashFreeDofs)
cbar = fig.colorbar(scalarMap)
#cbar = fig.colorbar(dots)
cbar.ax.set_yticks([-1,1])
cbar.update_ticks()
s=str(minDofs/float(maxDofs))
#cbar.set_label('mobility due to constraints')
cbar.set_label('remaining dofs')
#cbar.ax.set_yticklabels(['',''])
cbar.ax.set_yticklabels([str(minDofs),str(maxDofs)])
cbar.update_ticks()

plt.xlabel("RMSD to ini")
plt.ylabel("RMSD to target")
plt.savefig("pairedRMSD_twoTrees.png", transparent='True')
plt.clf()


os.chdir(currentPath)
d="energyPlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("energyPlots")

figEn = plt.figure()
ax = figEn.add_subplot(111, projection='3d')

ax.set_zlim(min(energy)-1000, max(energy)+1000)

ax.scatter(distanceToIni,distanceToTarget,energy, cmap=cmx.coolwarm, marker='s')

ax.set_xlabel('RMSD to ini')
ax.set_ylabel('RMSD to target')
ax.set_zlabel('Vdw energy')

plt.savefig("energy.png", transparent='True')
#plt.show()


orig_stdout = sys.stdout
fout = file('vdwEnergy.txt', 'w')
sys.stdout = fout

for i in range(len(energy)):
    print distanceToIni[i], distanceToTarget[i], energy[i]
    
sys.stdout = orig_stdout
fout.close()

if (len(reverse_energy) != 0):
    figEnRev = plt.figure()
    ax = figEnRev.add_subplot(111, projection='3d')
    
    ax.set_zlim(min(reverse_energy)-1000, max(reverse_energy)+1000)
    
    ax.scatter(reverse_distanceToTarget,reverse_distanceToIni,reverse_energy, cmap=cmx.coolwarm, marker='s')
    
    ax.set_xlabel('RMSD to ini')
    ax.set_ylabel('RMSD to target')
    ax.set_zlabel('Vdw energy')
    plt.savefig("reverse_energy.png", transparent='True')
    #plt.show()


    fout = file('vdwEnergyReverse.txt', 'w')
    sys.stdout = fout
    
    for i in range(len(reverse_energy)):
        print reverse_distanceToIni[i], reverse_distanceToTarget[i], reverse_energy[i]
        
    sys.stdout = orig_stdout
    fout.close()

os.chdir(currentPath)
