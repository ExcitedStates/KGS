#!/usr/bin/python
# coding: utf-8

# This file runs dilution analysis - rigidity analysis in KGS for a series of decreasing energy cut-off values.
# Start from the folder where all rigidity results shall be saved.

import numpy as np
import sys
import os
import click
import subprocess
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c,d,e): #,d,e,f
    # return a * np.exp(-b * x) + c #* np.exp(-d * x*x) +e
    return a*np.power(x,4)+b*np.power(x,3)+c*np.power(x,2)+d*x+e

def runningMean(x, N):
    return np.convolve(x, np.ones((N,))/N, mode='valid')

@click.command()
@click.option('--pdb', type=click.Path('rb'), required=True, help='PDB file')
@click.option('--initial', type=click.FLOAT, default= 0.0, help='Initial Mayo h-bond energy cutoff, default 0.0 kcal/mol.')
@click.option('--final', type=click.FLOAT, default= -6.1, help='Final Mayo h-bond energy cutoff, default -6.0 kcal/mol.')
@click.option('--step', type=click.FLOAT, default = 0.1, help='Mayo h-bond energy step size, default 0.1 kcal/mol.')
@click.option('--extrabonds', type=click.Path(exists=True), help='File with extra h-bonds to include as constraints (Format: ID1 ID2)')

def dilution(pdb, initial, final, step, extrabonds):
    
    currentDir = os.getcwd()
    #adapt to home location of ExcitedStates and repos
    baseFolder= currentDir[0:currentDir.find("ExcitedStates")] #if subfolder of the excited states repo
    
    rigidity= baseFolder+"KGS/build/kgs_rigidity --initial "
    baseName = pdb[:-4]
    print baseName
        
    hbondEnergies=np.arange(initial, final, -step)
    
    basicString = "python "+baseFolder+"KGS/scripts/kgs_prepare.py ../"+pdb
    # basicString = "kgs_prepare ../"+pdb
    # print basicString
    dihedrals=[]
    hbonds=[]
    cycleDofs=[]
    dofs=[]
    numClusters=[]
    biggestClusters=[]
    base=[]
    source=[]
    sink=[]
    shared=[]
    allostery=[]
    
    # # %%%%%%%%%%%% Computation of the whole thing
    # Run and process data
    for eCut in hbondEnergies:
        #Prepare file for KGS run
        # print "Preparing "+baseName+" with energy "+str(eCut)
        folderName = baseName+"_eCut"+str(eCut)
        if not os.path.exists(folderName):
            os.makedirs(folderName)
        os.chdir(folderName)
        kgsCall = basicString+" -energy "+str(eCut)+" -pre "+baseName+"_eCut"+str(eCut) # + " -noHydro"
        print kgsCall
        subprocess.call(kgsCall,shell=True)
        rigidityCall=rigidity+baseName+"_eCut"+str(eCut)+".kgs.pdb --workingDirectory ./ --source \"chain A\" --sink \"chain B\" > output.txt" #--source \"resi 180-191 and chain A\" --sink \"chain B and resi 149-169\" --hbondMethod user --hbondFile ../"+extrabonds+"
        #--source \"chain A\" --sink \"chain B\"
        subprocess.call(rigidityCall,shell=True)
        f = open("output/"+baseName+"_eCut"+str(eCut)+".kgs_stats_1.txt",'r')
        lines=f.readlines();
        count=0
        for line in lines:
            tokens = line.split(" ")
            if count==4:
                hbonds.append(int(tokens[4]))
            if count==8:
                dihedrals.append(int(tokens[6]))
            if count==10:
                cycleDofs.append(int(tokens[4]))
            elif count ==14:
                dofs.append(int(tokens[7]))
            elif count==19:
                numClusters.append(int(tokens[4]))
            elif count==20:
                biggestClusters.append(int(tokens[4]))
            count +=1
        f.close()
        f = open("output.txt",'r')
        lines=f.readlines();
        for line in lines:
            if line.startswith("Base"):
                tokens = line.split(" ")
                base.append(float(tokens[1]))
                source.append(float(tokens[3]))
                sink.append(float(tokens[5]))
                shared.append(float(tokens[7]))
                allostery.append( float(tokens[7])/(float(tokens[3])+float(tokens[5])) )
            if line.startswith("No coordinated motion"):
                base.append(0.0)
                source.append(0.0)
                sink.append(0.0)
                shared.append(0.0)
                allostery.append(0.0)
        f.close();
        os.chdir("../")
    with open('dilution.txt', 'w') as f:
        for item in range(0,len(numClusters)):
            f.write("%d %d %d %d %d %d %f %f %f %f %f\n" %(hbonds[item],dihedrals[item],cycleDofs[item],dofs[item],numClusters[item],biggestClusters[item],base[item],source[item],sink[item],shared[item],allostery[item]))
    
    # %%%%%%%%%%%% Computation of the whole thing
        
    # %%%%%%%%%%%% OLD READ IN <DEPRECATED FOR NEW SIMULATIONS
    # with open('dilution.txt', 'r') as f:
    #     for line in f.readlines():
    #         tokens=line.split(" ")
    #         hbonds.append(int(tokens[0]))
    #         dihedrals.append(int(tokens[1]))
    #         cycleDofs.append(int(tokens[2]))
    #         dofs.append(int(tokens[3]))
    #         numClusters.append(int(tokens[4]))
    #         biggestClusters.append(int(tokens[5]))
    #         source.append(float(tokens[6]))
    #         sink.append(float(tokens[7]))
    #         shared.append(float(tokens[8]))
    #         allostery.append(float(tokens[9]))
    # %%%%%%%%%%%% OLD READ IN <DEPRECATED FOR NEW SIMULATIONS
    
    # # # %%%%%%%%%%%% Read-in for quick plot when all data is already present
    # with open('dilution.txt', 'r') as f:
    #     for line in f.readlines():
    #         tokens=line.split(" ")
    #         hbonds.append(int(tokens[0]))
    #         dihedrals.append(int(tokens[1]))
    #         cycleDofs.append(int(tokens[2]))
    #         dofs.append(int(tokens[3]))
    #         numClusters.append(int(tokens[4]))
    #         biggestClusters.append(int(tokens[5]))
    #         base.append(float(tokens[6]))
    #         source.append(float(tokens[7]))
    #         sink.append(float(tokens[8]))
    #         shared.append(float(tokens[9]))
    #         allostery.append(float(tokens[10]))
    # # %%%%%%%%%%%%

    # Plot stuff
    fig, ax1 = plt.subplots()
    ax1.plot(hbondEnergies, hbonds,lw=4,label='h-bonds')
    ax1.plot(hbondEnergies, dihedrals,lw=4,label='|spanning tree|')
    ax1.plot(hbondEnergies, cycleDofs,lw=4,label='cycle DoF')
    ax1.plot(hbondEnergies, dofs,lw=4,label='internal DoF')
    ax1.plot(hbondEnergies, numClusters,lw=4,label='clusters')
    ax1.plot(hbondEnergies, biggestClusters,lw=4,label='|biggest cluster|')
    
    ax1.set_ylabel('count')
    ax1.set_xlabel('h-bond energy cutoff [kcal/mol]')
    # ax1.grid(True)
    
    # plot2Y = [x/(z1) for (x,z1) in zip(shared,dofs)]
    
    ax2 = ax1.twinx()
    # ax2.plot(hbondEnergies, allostery,lw=4,color='black',label='site DoF transfer')
    ax2.plot(hbondEnergies, allostery,lw=4,color='black',label='mutual information')
    ax2.set_ylabel('mutual information [bits]')
    ax1.legend(loc=0)#,bbox_to_anchor=(0.0, 0.5, 0.5, 0.5))
    ax2.legend(loc=6,bbox_to_anchor=(0.0, 0., 0.5, 0.5))
    plt.savefig("dilutionPlot.eps")
    plt.savefig("dilutionPlot.png")
    #plt.show()
    plt.clf()
    
    ## Floppy mode density and analysis
    floppyModeDensity = np.asarray([float(x)/float(y) for x,y in zip(dofs,dihedrals)])
    popt, pcov = curve_fit(func, hbondEnergies, floppyModeDensity)    
    
    #Derivatives from moving mean
    # misFit = 3;
    # dataFit = runningMean(floppyModeDensity,3)
    #Derivatives from fitted curves
    # dataFit = func(hbondEnergies,*popt)
    
    #Derivatives from original data or fitted curve
    firstDerivFMD = np.gradient(floppyModeDensity) # 2nd order central diff scheme on interior points, 1st order on boundary points
    # firstDerivFMD = np.gradient(dataFit,edge_order=2) # 2nd order central diff scheme on interior points, 1st order on boundary points
    
    secondDerivFMD = np.gradient(firstDerivFMD,edge_order=2) # 2nd order central diff scheme on interior points, 1st order on boundary points
    
    fig, ax1 = plt.subplots()
    ax1.plot(hbondEnergies, floppyModeDensity,lw=4,color='royalblue')
    # ax1.plot(hbondEnergies[1:-1], dataFit,lw=4,color='aqua')
    
    ax1.set_xlabel('h-bond energy cutoff [kcal/mol]')
    ax1.set_ylabel('floppy mode density (fmd)',color='royalblue')
    ax1.tick_params('y', colors='royalblue')
    
    ax2 = ax1.twinx()
    ax2.plot(hbondEnergies, secondDerivFMD,lw=4,color='lightcoral')
    ax2.set_ylabel('second derivative', color='lightcoral')
    ax2.tick_params('y', colors='lightcoral')
    
    # plt.grid(True)
    fig.tight_layout()
    plt.savefig("dilutionPlotPhi.eps")
    plt.savefig("dilutionPlotPhi.png")
    #plt.show()
    plt.clf()
    
        
if __name__ == "__main__":
    dilution()