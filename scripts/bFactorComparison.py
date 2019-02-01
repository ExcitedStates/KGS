#!/usr/bin/env python

import sys
import math
import numpy as np
import numpy.random
from sklearn.decomposition import PCA
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from pdb_structure import *
import pylab
import scipy.stats as st


def consistentResidueList(pdbFileNames):
	#Returns a list of residue IDs that are present in all structure files (checked pairs of name/ID)
	residueList=[]
	first=True
	for fn in pdbFileNames:
		pdbFile = PDBFile(fn)
		thisList = pdbFile.getResidueIDsandNames() #List of pairs of residue IDs and residue names
		if first:
			first=False
			residueList = thisList
		residueList = [i for i in thisList if i in residueList]
	residueList = [entry[0] for entry in residueList]
	return residueList

def bFactors(pdbFileNames, residueList):
	X=[]
	for fn in pdbFileNames:
		pdbFile = PDBFile(fn)
		#Pairwise atomic distance matrix for each structure
		bFacs = pdbFile.bFactorList(0,["CA"],residueList)
		# printVecToFile(atomMat,fn[:-4]+"_bFac.txt")
		n = len(bFacs)
		#Convert atomic-distance matrix to observation vector
		# observationVector=[]
		# for i in range(n-1):
		# 	for j in range(i+1, n):
		# 		observationVector.append(atomMat[i][j])
		#Save observation as one row in X
		# X.append(observationVector)
		X.append(bFacs)
	return X

if __name__ == "__main__":
	if(len(sys.argv)<=2):
		print "Usage:",sys.argv[0],"<pdb-files in a row>"
		sys.exit(0)
		
	pdbFileNames = sys.argv[1:]
	# self.pdbFileNames = [ x[x.rfind("/")+1:x.rfind(".pdb")] for x in pdbFileNames]
	# self.pdbFileNames = [ x.replace(".kgs","") for x in self.pdbFileNames]
	indices = range(len(pdbFileNames))
	residueList = consistentResidueList(pdbFileNames)
	# labels=["conv_A","conv_B","xfel_A","xfel_B","raxis_A","raxis_B"]
	# labels=["structure 1","structure 2"]
	labels=[fileName[0:-4] for fileName in pdbFileNames]
	
	X = bFactors(pdbFileNames,residueList)
	# f=plt.figure()
	for ent1 in range(len(X)):
		fn = pdbFileNames[ent1]
		plt.plot(X[ent1],label=labels[ent1],alpha=0.7,lw=2)
		print "Mean:",np.mean(X[ent1])
		# cbar=plt.colorbar()
	plt.xlabel('residue ID')
	plt.ylabel('B factor')
	plt.legend()
	plt.savefig("bFactorComparison.png", transparent='False',format='png',dpi=600)
	plt.close()
	
	for ent1 in range(len(X)/2):
		studTVal = st.ttest_ind(X[ent1*2],X[ent1*2+1],equal_var=True)
		# welchTVal = st.ttest_ind(X[ent1*2],X[ent1*2+1],equal_var=False)
		# pVal = st.pearsonr(X[ent1*2],X[ent1*2+1])
		# print "Pearson correlation ",pVal[0],pVal[1]
		print "Independent Student's t-test: ",studTVal
		# print "Non-equal variance Welch's t-test: ",welchTVal