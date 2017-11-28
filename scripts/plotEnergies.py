#!/usr/bin/python

import sys
import os
import numpy as np
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
from matplotlib.ticker import MultipleLocator, LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from mpl_toolkits.mplot3d import Axes3D


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from pylab import rcParams
rcParams['figure.figsize'] = 6.5, 4

ftSize=12
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

if len(sys.argv)==1:
    print "Usage: plotEnergies.py ", sys.argv[0]," <output directory>"

enablelog = True

rmsdIni=[]
rmsdTarget=[]
energy=[]
# reverseRmsdIni=[]
# reverseRmsdTarget=[]
# reverseEnergy=[]
clashFreeDofsForward=[]
clashFreeDofsReverse=[]

relevant_path = sys.argv[1]
included_extensions = ['*[0-9].pdb'] ;

all_files =[fn for fn in os.listdir(relevant_path) if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) ];
file_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and not fn.startswith('_')];
target_names = [fn for fn in all_files if any ([fnmatch.fnmatch(fn,ext) for ext in included_extensions]) and fn.startswith('_')];

regPlanner=1
maxDofs=0
minDofs=99999
dofs=0

fwdYLimit=9999
fwdXLimit=0.0
fwdMaxDist=0.0
revYLimit=0
revXLimit=99999
revMaxDist=0.0
iniDistance=1.0

with open("output.txt") as outputFile:
    for line in outputFile:
        if "Mean distance to target," in line:
            tokens = line.split(' ')
            iniDistance = float(tokens[11])
            fwdMaxDist = iniDistance
            revMaxDist = iniDistance

        if "Using regular planner" in line:
            regPlanner = 1
        if "Switching to reversePlanner" in line:
            regPlanner = 0
        if "RMSD to target:" in line:
            tokens = line.split(' ')
            rmsdTarget =tokens[3]
            rmsdTarget = float(rmsdTarget[:-1])
            if(regPlanner == 1):
                if(rmsdTarget < fwdYLimit):
                    fwdYLimit = rmsdTarget
                if(rmsdTarget > fwdMaxDist):
                    fwdMaxDist = rmsdTarget
            if (regPlanner == 0):
                if(rmsdTarget < revXLimit):
                    revXLimit = rmsdTarget
                if(rmsdTarget > revMaxDist):
                    revMaxDist = rmsdTarget
        if "New structure" in line:
            tokens = line.split(' ')
            dofs =int(tokens[21])
            if(regPlanner == 1):
                clashFreeDofsForward.append(dofs)
            else:
                clashFreeDofsReverse.append(dofs)
            if( dofs > maxDofs):
                maxDofs = dofs
            if( dofs < minDofs):
                minDofs=dofs        
                
currentPath = os.getcwd()
os.chdir(relevant_path)

for pdb_file_name in file_names:
        rmsdToIni = 0.0
        rmsdToTarget = 1000.0
        en=0.0
        with open(pdb_file_name) as pdb_file:
            for line in pdb_file:
                if "REMARK\tRMSD to initial" in line:
                    rmsdToIni = float(line[25:])
                if "REMARK\tRMSD to target" in line:
                    rmsdToTarget = float(line[24:])
                if "REMARK\tVdw energy" in line:
                    en = float(line[20:])

        rmsdIni.append(rmsdToIni)
        rmsdTarget.append(rmsdToTarget)
        energy.append(en)

switch=len(rmsdIni)
        
for pdb_file_name in target_names:
        reverseRmsdToIni = 0.0
        reverseRmsdToTarget = 1000.0
        with open(pdb_file_name) as pdb_file:
            for line in pdb_file:
                if "REMARK\tRMSD to initial" in line:
                    reverseRmsdToIni = float(line[25:])
                if "REMARK\tRMSD to target" in line:
                    reverseRmsdToTarget = float(line[24:])
                if "REMARK\tVdw energy" in line:
                    en = float(line[20:])

        # reverseRmsdIni.append(reverseRmsdToIni)
        # reverseRmsdTarget.append(reverseRmsdToTarget)
        # reverseEnergy.append(en)

        rmsdIni.append(reverseRmsdToTarget)
        rmsdTarget.append(reverseRmsdToIni)
        energy.append(en)

os.chdir(currentPath)
d="energyPlots"
if not os.path.exists(d):
    os.makedirs(d)
os.chdir("energyPlots")

fig = plt.figure()

# ---------------------------
#3D Plot

# ax = fig.add_subplot(111, projection='3d')
# # ax.set_zscale('log')
# # surf = ax.plot_surface(rmsdIni,rmsdTarget,energy, rstride=1, cstride=1, cmap=cm.coolwarm,
# #         linewidth=0, antialiased=False)

# ax.set_zlim(min(energy)-1000, max(energy)+1000)
# # ax.zaxis.set_major_locator(LinearLocator(10))
# # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# # fig.colorbar(surf, shrink=0.5, aspect=5)

# ax.scatter(rmsdIni,rmsdTarget,energy, cmap=cm.coolwarm, marker='s')

# ax.set_xlabel('RMSD to ini')
# ax.set_ylabel('RMSD to target')
# ax.set_zlabel('Vdw energy')

# # plt.xlabel("RMSD to ini")
# # plt.ylabel("RMSD to target")
# plt.savefig("energy.png")
# plt.show()
# ---------------------------

# ---------------------------
#@D Plot with color

ax = fig.add_axes()

minVal = min(energy)
maxVal = max(energy)
jet = cm = plt.get_cmap('jet')
#jet = cm = plt.get_cmap('coolwarm') 

cNorm  = colors.Normalize(vmin=minVal, vmax=maxVal)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

for i in range(len(energy)):
    colorval=scalarMap.to_rgba(energy[i])
    #colorval=scalarMap.to_rgba(clashFreeDofsForward[i])
    if i<=switch:
        dot, = ax.plot(rmsdIni[i],rmsdTarget[i],'v',color=colorval,markeredgecolor=colorval,alpha=0.7)
    else:
        dot, = ax.plot(rmsdIni[i],rmsdTarget[i],'^',color=colorval,markeredgecolor=colorval,alpha=0.7)
    dots.append(dot)

ax.set_xlabel('RMSD to ini')
ax.set_ylabel('RMSD to target')

scalarMap.set_array(dots)

plt.savefig("energyLandscape.png")


orig_stdout = sys.stdout
fout = file('vdwEnergy.txt', 'w')
sys.stdout = fout

for i in range(len(energy)):
    print rmsdIni[i], rmsdTarget[i], energy[i]
    
sys.stdout = orig_stdout
fout.close()

fout = file('vdwEnergyReverse.txt', 'w')
sys.stdout = fout

for i in range(len(reverseEnergy)):
    print reverseRmsdIni[i], reverseRmsdTarget[i], reverseEnergy[i]
    
sys.stdout = orig_stdout
fout.close()

os.chdir(currentPath)
