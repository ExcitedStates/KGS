#!/usr/bin/python

import sys
import os
import numpy as np
from combineKGSPath_steps import extractPath
from operator import add

def main():

    if len(sys.argv)<2:
        print "Usage: "+sys.argv[0]+" <path.pdb files in a row> "
        sys.exit(1)

    dof_ids=[]
    resids=[]
    dof_kinds=[]
    torsions=[]
    sumSteps=[]
    mainChain=[]
    torsion=[]
    sumStep=[]

    for i in range(len(sys.argv)-1):
        pdbPath=sys.argv[i+1]
        pathList, reversePathList = extractPath(pdbPath)
        finalConf=pathList[-1]
        print "Final conf: "+str(finalConf)
        pathFileSepIdx = pdbPath.rfind("/")
        outputDir = pdbPath[0:pathFileSepIdx] if pathFileSepIdx!=-1 else "."
        # sampleNum = int(pdbPath[pdbPath.rfind("_")+1:pdbPath.rfind(".pdb")])
        modelName = str(pdbPath[pdbPath.rfind("/")+1:pdbPath.rfind("_path")])
        
        currentTorsions=[]
        currentSum=[]
        # fileString=outputDir+"/"+modelName+"_q_"+str(finalConf)+"_qCombo.txt"
        fileString = outputDir+"/"+modelName+"_q_"+str(finalConf)+".txt"
        print "opening "+fileString
        with open(fileString, "r") as q_file:
            for line in q_file:
                tokens = line.split(' ')

                dof_ids.append(int(tokens[0]))
                resids.append(int(tokens[1]))
                dof_kinds.append(int(tokens[2]))
                torsions.append(float(tokens[3]))
                sumSteps.append(float(tokens[4]))
                mainChain.append(int(tokens[5]))

        # for j in pathList:
        # 	fileString = outputDir+"/"+modelName+"_q_"+str(j)+".txt"
        # 	with open(fileString, "r") as q_file:
	       #      for line in q_file:
	       #          tokens = line.split(' ')
	       #          torsion.append(float(tokens[3]))
	       #          sumStep.append(float(tokens[4]))

	       #  for k in range(len(torsion)):
	       #  	torsions[k] = torsions[k]+abs(torsion[k])
    
    d="comboAnalysis"
    if not os.path.exists(d):
        os.makedirs(d)
    os.chdir("comboAnalysis")
    
    orig_stdout = sys.stdout
    fout = file('combined_qFile.txt', 'w')
    sys.stdout = fout

    for i in range(len(dof_ids)):
        print dof_ids[i], resids[i], dof_kinds[i], torsions[i], sumSteps[i], mainChain[i]

    sys.stdout = orig_stdout
    fout.close()
    
if __name__ == "__main__":
	main()