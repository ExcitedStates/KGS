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
# from bFactorComparison import bFactors

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

def atomDistanceMatrix(pdbFile,residueList):
	coordMat = pdbFile.coordMatrix(0,["CA"],residueList) #adapt to use with other atom types as well
	atoms = coordMat.shape[0]
	assert atoms>0

	dmat = np.zeros(  shape=( atoms,atoms )  )
	#Matrix of distance norms between atoms in one Structure
	for i in range(atoms-1):
		for j in range(i+1,atoms):
			d = coordMat[i]-coordMat[j]
			dmat[i][j] = np.linalg.norm(d)

	return dmat

def differenceDataMatrix(pdbFileNames, residueList):
	X=[]
	for fn in pdbFileNames:
		pdbFile = PDBFile(fn)
		#Pairwise atomic distance matrix for each structure
		atomMat = atomDistanceMatrix(pdbFile,residueList)
		printMatToFile(atomMat,fn[:-4]+"_dMat.txt")
		n = len(atomMat)
		#Convert atomic-distance matrix to observation vector
		# observationVector=[]
		# for i in range(n-1):
		# 	for j in range(i+1, n):
		# 		observationVector.append(atomMat[i][j])
		#Save observation as one row in X
		# X.append(observationVector)
		X.append(atomMat)
	diffMat=[]
	# for ent in range(int(0.5*len(X[0])*(len(X[0])-1))):
	for ent1 in range(len(X)-1):
		for ent2 in range(ent1+1,len(X)):
			fn1=pdbFileNames[ent1]
			fn2=pdbFileNames[ent2]
			dMat = X[ent1]-X[ent2]
			printMatToFile(dMat,fn1[:-4]+"_"+fn2[:-4]+"_diffMat.txt")
			diffMat.append(dMat)
	return diffMat

def printMatToFile(X,fileName):
	orig_stdout = sys.stdout
	fout = file(fileName,'w')
	sys.stdout = fout
	for row in X:
		print '\t'.join(map(str,row))
	sys.stdout = orig_stdout
	fout.close()

	
if __name__ == "__main__":
	if(len(sys.argv)<=2):
		print "Usage:",sys.argv[0],"<pdb-file 1> <pdb-file 2>"
		sys.exit(0)
		
	pdbFileNames = sys.argv[1:]
	# self.pdbFileNames = [ x[x.rfind("/")+1:x.rfind(".pdb")] for x in pdbFileNames]
	# self.pdbFileNames = [ x.replace(".kgs","") for x in self.pdbFileNames]
	indices = range(len(pdbFileNames))
	residueList = consistentResidueList(pdbFileNames)
	X = differenceDataMatrix(pdbFileNames,residueList)
	for ent1 in range(len(X)):
		fn = pdbFileNames[ent1]
		plt.matshow(np.transpose(np.array(X[ent1])),cmap='bwr')
		cbar=plt.colorbar()
		cbar.set_label('difference distance [Angstrom]')
		plt.savefig("diffMat"+fn[:-4]+"_"+str(ent1)+".png", transparent='False',format='png',dpi=600)
		plt.close()
		
	# X = bFactors(pdbFileNames,residueList)
	# for ent1 in range(len(X)):
	# 	fn = pdbFileNames[ent1]
	# 	plt.plot(X[ent1])
	# 	print "Mean:",np.mean(X[ent1])
	# 	# cbar=plt.colorbar()
	# # plt.set_label('B factor [Angstrom]')
	# plt.savefig("bFacComp.png", transparent='False',format='png',dpi=600)
	# plt.close()
		