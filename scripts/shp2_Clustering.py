#!/usr/bin/python

import sys
import math
import numpy as np
import numpy.random
import scipy.cluster.hierarchy
from pdb_structure import *


class KMedoidsClustering:
	def __init__(self, pdbFileNames, k=10):
		self.pdbFileNames = [ x.replace(".pdb","") for x in pdbFileNames]
		n = len(pdbFileNames)

		#Initial clustering based on atom-based RMSD in structures
		RMSDmatrix = HierarchicalClustering.dRMSDmatrix(pdbFileNames)
		
		# initialize k random cluster centers (mediods) from the n structures
		self.medoids = list(numpy.random.randint(0,n,k)) #indicated by indices
		#ToDo: check what happens if we start from hierarchical centroids for the clusters
		#Compute cost (distance from each sample to closest medoid)
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
					# newCost = KMedoidsClustering.medoidCost(self.medoids, subsample, RMSDmatrix)
					newCost = KMedoidsClustering.medoidCost(self.medoids, range(n), RMSDmatrix)
					if newCost<cost:
						print 'Improved cost:',cost,'->',newCost
						cost = newCost
						break
					else:
						self.medoids.remove(o)
						self.medoids.append(m)

			cost = KMedoidsClustering.medoidCost(self.medoids, range(n), RMSDmatrix)
			print 'Real cost',cost
			
		clusterList=[]
		for i in range(n):
			row = RMSDmatrix[i][:]
			center = np.argmin([row[m] for m in self.medoids])
			clusterList.append((i,self.medoids[center]))
		self.clusterList = clusterList


	@staticmethod
	def medoidCost(medoids, subsample, dMat):
		cost = 0.0
		for i in subsample:
			cost += min([dMat[i][m] for m in medoids]) #only picks cost to closest medoid
		return cost/len(subsample)
		#This also contains the medoids themselves in the samples, but their min distance is zero
		
	
class HierarchicalClustering:

	def __init__(self, pdbFileNames, desired_clusters=10):
		"""
		Separates the PDB-files into the desired number of clusters based on dRMSD differences. The clustering
		is stored in the self.clusterIndices field which contains len(pdbFileNames) entries indicating the cluster-id
		of each pdb file (1-indexed).
		"""
		self.pdbFileNames = pdbFileNames

		#RMSD across structures computed from inter-atomic distance matrices within each structure
		RMSDmatrix = HierarchicalClustering.dRMSDmatrix(pdbFileNames)

		n = len(pdbFileNames)
		#Convert RMSDmatrix to array containing only the upper-triangle values (input format for scipy HC)
		RMSDvector = []
		for i in range(n-1):
			for j in range(i+1, n):
				RMSDvector.append(RMSDmatrix[i][j])

		print 'Clustering .. '
		#Hierarchical clustering based on RMSD matrix
		cl = scipy.cluster.hierarchy.complete(RMSDvector)
		self.clusterIndices = scipy.cluster.hierarchy.fcluster(cl, desired_clusters, criterion='maxclust')
		print 'done'
		# print fcl
		#Compute the root nodes (L) in the hierarchical clustering
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
	def atomDistanceMatrix(pdbFile,residueList):
		coordMat = pdbFile.coordMatrix(0,["CA"],residueList)
		atoms = coordMat.shape[0]
		assert atoms>0

		dmat = np.zeros(  shape=( atoms,atoms )  )
		#Matrix of distance norms between atoms in one Structure
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
		
		#Preprocessing step to identify residues present in each pdb
		print 'Preprocessing residue information'
		residueList=[]
		first=True
		for fn in pdbFileNames:
			# if i%(max(100,len(pdbFileNames))/100)==0:
			# 	sys.stdout.write('.')
			# 	sys.stdout.flush()
			# i+=1

			pdbFile = PDBFile(fn)
			thisList = pdbFile.getResidueIDsandNames() #List of pairs of residue IDs and residue names
			if first:
				first=False
				residueList = thisList
			residueList = [i for i in thisList if i in residueList]
		residueList = [entry[0] for entry in residueList]	
		
		print 'Reading files'
		i = 0
		for fn in pdbFileNames:
			if i%(max(100,len(pdbFileNames))/100)==0:
				sys.stdout.write('.')
				sys.stdout.flush()
			i+=1

			pdbFile = PDBFile(fn)
			#Pairwise atomic distance matrix for each structure
			pdbtoamat_map[ fn ] = HierarchicalClustering.atomDistanceMatrix(pdbFile,residueList)
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
				#RMSD of pairwise atomic distances across structures
				rmsd = HierarchicalClustering.dRMSD( amat_i, amat_j )
				matrix[i][j] = rmsd
				matrix[j][i] = rmsd
		print ' done'
		return matrix

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
				print "Computing "+str(clusters)+" clusters"
			except:
				print "Commandline clusters format: -clusters <int>"
				sys.exit(-1)
			break


	if(len(sys.argv)<=1):
		print "Usage:",sys.argv[0],"[-clusters <int>] [-verbose] <list of pdb-files>"
		sys.exit(0)
	# hc = HierarchicalClustering(sys.argv[1:])
	kc = KMedoidsClustering(sys.argv[1:],clusters)

	# if verbose:
	# 	print hc.clusterIndices
	print kc.clusterList
	count=1
	for c in kc.medoids:
		memberList=[kc.pdbFileNames[x[0]] for x in kc.clusterList if x[1] == c ]
		# print "Center "+kc.pdbFileNames[c]
		# print "Members "+' '.join(memberList)
		print "select cluster"+str(count)+", "+'* or '.join(memberList)
		count += 1
		for mem in memberList:
			print "align "+mem+"*, "+kc.pdbFileNames[c]+"*"
			# print "align "+mem+"* and PTP, "+kc.pdbFileNames[c]+"* and PTP"