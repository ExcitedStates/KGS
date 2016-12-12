#!/usr/bin/python

import sys
import math
import numpy as np
import numpy.random
import scipy.cluster.hierarchy
from pdb_structure import *


class KMedoidsClustering:
	def __init__(self, pdbFileNames, k=10):
		self.pdbFileNames = pdbFileNames
		n = len(pdbFileNames)

		RMSDmatrix = HierarchicalClustering.dRMSDmatrix(pdbFileNames)

		self.medoids = list(numpy.random.randint(0,n,k))

		cost = KMedoidsClustering.medoidCost(self.medoids, range(n), RMSDmatrix)
		cost_convergence = 1000000
		while cost<cost_convergence:
			cost_convergence = cost
			for m in self.medoids:
				for o in range(n):
#				subsample = numpy.random.randint(0,n,400) if n>400 else range(n)
#				cost = KMedoidsClustering.medoidCost(self.medoids, subsample, RMSDmatrix)
#				for o in subsample:
					self.medoids.remove(m)
					self.medoids.append(o)
					newCost = KMedoidsClustering.medoidCost(self.medoids, subsample, RMSDmatrix)
					#newCost = KMedoidsClustering.medoidCost(self.medoids, range(n), RMSDmatrix)
					if newCost<cost:
						print 'Improved cost:',cost,'->',newCost
						cost = newCost
						break
					else:
						self.medoids.remove(o)
						self.medoids.append(m)

			cost = KMedoidsClustering.medoidCost(self.medoids, range(n), RMSDmatrix)
			print 'Real cost',cost



	@staticmethod
	def medoidCost(medoids, subsample, dMat):
		cost = 0.0
		for i in subsample:
			cost += min([dMat[i][m] for m in medoids])
		return cost/len(subsample)


class HierarchicalClustering:

	def __init__(self, pdbFileNames, desired_clusters=10):
		"""
		Separates the PDB-files into the desired number of clusters based on dRMSD differences. The clustering
		is stored in the self.clusterIndices field which contains len(pdbFileNames) entries indicating the cluster-id
		of each pdb file (1-indexed).
		"""
		self.pdbFileNames = pdbFileNames

		RMSDmatrix = HierarchicalClustering.dRMSDmatrix(pdbFileNames)

		n = len(pdbFileNames)
		#Convert RMSDmatrix to array containing only the upper-triangle values (input format for scipy HC)
		RMSDvector = []
		for i in range(n-1):
			for j in range(i+1, n):
				RMSDvector.append(RMSDmatrix[i][j])

		print 'Clustering .. '
		cl = scipy.cluster.hierarchy.complete(RMSDvector)
		self.clusterIndices = scipy.cluster.hierarchy.fcluster(cl, desired_clusters, criterion='maxclust')
		print 'done'
		# print fcl
		L,M = scipy.cluster.hierarchy.leaders(cl,self.clusterIndices)

		self.centroids = []
		for centroidIdx in L:
			self.centroids.append( pdbFileNames[centroidIdx] )
		# print L
		# print M


		assert len(self.clusterIndices)==len(pdbFileNames)
		assert min(self.clusterIndices)>0
		assert max(self.clusterIndices)<=desired_clusters


	@staticmethod
	def atomDistanceMatrix(pdbFile):
		coordMat = pdbFile.coordMatrix(names=["C1'","CB"])
		atoms = coordMat.shape[0]
		assert atoms>0

		dmat = np.zeros(  shape=( atoms,atoms )  )
		for i in range(atoms-1):
			for j in range(i+1,atoms):
				d = coordMat[i]-coordMat[j]
				dmat[i][j] = np.linalg.norm(d)
			# dmat[j][i] = dmat[i][j]

		return dmat

	@staticmethod
	def dRMSD(mat1, mat2):
		assert mat1.shape == mat2.shape
		assert mat1.shape[0]==mat1.shape[1]

		atoms = mat1.shape[0]
		mat3 = np.square(mat2-mat1)
		return math.sqrt(mat3.sum()/float(atoms*(atoms-1)))


	@staticmethod
	def dRMSDmatrix(pdbFileNames):
		n = len(pdbFileNames)
		assert n>0

		matrix = np.zeros(shape=(n,n))

		pdbtoamat_map = {} #Maps PDB file name to atom distance matrix
		print 'Reading files'
		i = 0
		for fn in pdbFileNames:
			if i%(max(100,len(pdbFileNames))/100)==0:
				sys.stdout.write('.')
				sys.stdout.flush()
			i+=1

			pdbFile = PDBFile(fn)
			pdbtoamat_map[ fn ] = HierarchicalClustering.atomDistanceMatrix(pdbFile)
		print ' done'
		print 'Filling matrix'
		#Fill matrix
		for i in range(n-1):
			if i%(max(n,100)/100)==0:
				sys.stdout.write('.')
				sys.stdout.flush()

			# print '>',i,'of',n
			amat_i = pdbtoamat_map[pdbFileNames[i]]
			for j in range(i+1,n):
				amat_j = pdbtoamat_map[pdbFileNames[j]]

				rmsd = HierarchicalClustering.dRMSD( amat_i, amat_j )
				matrix[i][j] = rmsd
				matrix[j][i] = rmsd
		print ' done'
		return matrix




# pdbFiles = ['/Users/rfonseca/Documents/KGSDock/BenchmarkSet/Experiments/2B6G/rna_sampling/ri20_ns_r90_a2/output/newpdb_'+str(i)+'.pdb' for i in range(1,50)]
# km = KMedoidsClustering(pdbFiles)
# sys.exit(0)
#
# hc = HierarchicalClustering(pdbFiles)
# amat1 = HierarchicalClustering.atomDistanceMatrix(PDBFile(pdbFiles[0]))
# amat2 = HierarchicalClustering.atomDistanceMatrix(PDBFile(pdbFiles[1]))
# print HierarchicalClustering.dRMSD(amat1, amat2)

if __name__ == "__main__":

	verbose = False
	for i in range(len(sys.argv)):
		if sys.argv[i]=='-verbose':
			del sys.argv[i]
			verbose = True
			break


	clusters = 10
	for i in range(len(sys.argv)):
		if sys.argv[i]=='-clusters':
			try:
				clusters = int(sys.argv[i+1])
				del sys.argv[i:i+2]
			except:
				print "Commandline clusters format: -clusters <int>"
				sys.exit(-1)
			break


	if(len(sys.argv)<=1):
		print "Usage:",sys.argv[0],"[-clusters <int>] [-verbose] <list of pdb-files>"
		sys.exit(0)
	# hc = HierarchicalClustering(sys.argv[1:])
	kc = KMedoidsClustering(sys.argv[1:])

	# if verbose:
	# 	print hc.clusterIndices

	for c in kc.medoids:
		print kc.pdbFileNames[c]
