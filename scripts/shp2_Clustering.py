#!/usr/bin/python

import sys
import math
import numpy as np
import numpy.random
import scipy.cluster.hierarchy
from sklearn import cluster, datasets
from sklearn.decomposition import PCA
from pdb_structure import *
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt

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

	def clusterSilhouette(self):
		for cluster in self.clusterList:
			for member in cluster:
				pass
	
	def printClusterInfo(self):
		print self.clusterList
		colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
		'salmon', 'smudge', 'hotpink', 'slate', 'magenta', 'orange', 'marine', \
		'olive', 'density', 'teal', 'forest', 'firebrick', 'chocolate',    \
		'wheat', 'white']
		cluster=1
		colorlist=[]
		for c in self.medoids:
			memberList=[self.pdbFileNames[x[0]] for x in self.clusterList if x[1] == c ]
			print "select cluster"+str(cluster)+", "+'* or '.join(memberList)
			print "color "+colours[cluster-1]+", cluster"+str(cluster)
			cluster += 1
			for mem in memberList:
				print "align "+mem+"*, "+self.pdbFileNames[c]+"*"
				# print "align "+mem+"* and PTP, "+self.pdbFileNames[c]+"* and PTP"
	
class HierarchicalClustering:

	def __init__(self, pdbFileNames, desired_clusters=10):
		"""
		Separates the PDB-files into the desired number of clusters based on dRMSD differences. The clustering
		is stored in the self.clusterIndices field which contains len(pdbFileNames) entries indicating the cluster-id
		of each pdb file (1-indexed).
		"""
		self.pdbFileNames = [ x.replace(".pdb","") for x in pdbFileNames]

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
		cl = scipy.cluster.hierarchy.linkage(RMSDvector,method='complete') #complete --> furthest distance between clusters; 'single' --> shortest distance
		self.clusterIndices = scipy.cluster.hierarchy.fcluster(cl, desired_clusters, criterion='maxclust') #maximum of desired_clusters
		print 'done'
		print self.clusterIndices

		assert len(self.clusterIndices)==len(pdbFileNames)
		assert min(self.clusterIndices)>0
		assert max(self.clusterIndices)<=desired_clusters
		
		#Compute medoids for the clusters
		medoids = HierarchicalClustering.computeMedoids(self.clusterIndices,desired_clusters,RMSDmatrix)
		self.medoids=[]
		for medoidID in medoids:
			self.medoids.append(self.pdbFileNames[medoidID])
		print self.medoids
		
	def printClusterInfo(self):
		colours = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
		'salmon', 'smudge', 'hotpink', 'slate', 'magenta', 'orange', 'marine', \
		'olive', 'density', 'teal', 'forest', 'firebrick', 'chocolate',    \
		'wheat', 'white']
		cluster=1
		for c in self.medoids:
			memberList=[self.pdbFileNames[x[0]] for x in enumerate(self.clusterIndices) if x[1]==cluster ]
			print "select cluster"+str(cluster)+", "+'* or '.join(memberList)
			print "color "+colours[cluster-1]+", cluster"+str(cluster)
			cluster += 1
			for mem in memberList:
				print "align "+mem+"*, "+c+"*"
				# print "align "+mem+"* and PTP, "+self.pdbFileNames[c]+"* and PTP"
				
	@staticmethod
	def computeMedoids(clusterIndices,desired_clusters,RMSDmatrix):
		medoids=[]
		print "Identifying best medoids in hierarchical clusters"
		for cluster in range(desired_clusters):
			members = [x[0] for x in enumerate(clusterIndices) if x[1]==cluster+1 ]
			cost = 1000000
			bestMedoid=0
			for medoid in members:
				newCost = KMedoidsClustering.medoidCost([medoid],members,RMSDmatrix)
				if newCost<cost:
					print 'Improved cost:',cost,'->',newCost
					cost = newCost
					bestMedoid = medoid
			medoids.append(bestMedoid)
		return medoids

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
		residueList=consistentResidueList(pdbFileNames)	
		
		print 'Reading files'
		i = 0
		for fn in pdbFileNames:
			if i%(max(100,len(pdbFileNames))/100)==0:
				sys.stdout.write('.')
				sys.stdout.flush()
			i+=1

			pdbFile = PDBFile(fn)
			#Pairwise atomic distance matrix for each structure
			pdbtoamat_map[ fn ] = atomDistanceMatrix(pdbFile,residueList)
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

def torsion(a1, a2, a3, a4):
	""" 
	Taken from 
	http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
	"""
	from numpy import cross, dot
	from numpy.linalg import norm
	import math
	b1 = (a2-a1)
	b2 = (a3-a2)
	b3 = (a4-a3)
	n1 = cross(b1,b2)/norm( cross(b1,b2) )
	n2 = cross(b2,b3)/norm( cross(b2,b3) )
	m1 = cross( n1, b2/norm(b2) )
	x = dot(n1, n2)
	y = dot(m1, n2)
	angle = -math.atan2(y,x)
	return angle

def backboneTorsionMatrix(pdbFileNames,residueList):
	X=[]
	for fn in pdbFileNames:
		pdbFile = PDBFile(fn)
		#Backbone torsions in vector form
		X.append(backboneTorsionVector(pdbFile,residueList))
	return X

def backboneTorsionVector(pdbFile,residueList):
	torsions=[]
	triplets = [residueList[i:i+3] for i in range(len(residueList)-2)]
	for entry in triplets:
		if ((entry[0]+1 != entry[1]) or (entry[1]+1 != entry[2])):
			continue
		a1 = pdbFile.getAtom(entry[0] , "C")
		a2 = pdbFile.getAtom(entry[1] , "N")
		a3 = pdbFile.getAtom(entry[1] , "CA")
		a4 = pdbFile.getAtom(entry[1] , "C")
		a5 = pdbFile.getAtom(entry[2] , "N")
		if not (a2 and a3 and a4):
			continue;
		if not a1:
			continue;
		phi = torsion(a1,a2,a3,a4)
		torsions.append(phi)
		if not a5:
			continue
		psi  = torsion(a2,a3,a4,a5)
		torsions.append(psi)
	return torsions
	
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

	return dmat
	
def fullDataMatrix(pdbFileNames, residueList):
	X=[]
	for fn in pdbFileNames:
		pdbFile = PDBFile(fn)
		#Pairwise atomic distance matrix for each structure
		atomMat = atomDistanceMatrix(pdbFile,residueList)
		n = len(atomMat)
		#Convert atomic-distance matrix to observation vector
		observationVector=[]
		for i in range(n-1):
			for j in range(i+1, n):
				observationVector.append(atomMat[i][j])
		#Save observation as one row in X
		X.append(observationVector)
	return X

def printMatToFile(X,fileName):
	orig_stdout = sys.stdout
	fout = file(fileName,'w')
	sys.stdout = fout
	for row in X:
		print '\t'.join(map(str,row))
	sys.stdout = orig_stdout
	fout.close()
	
class DBscan:

	def __init__(self, pdbFileNames):
		"""
		Separates the PDB-files into the desired number of clusters based on interatomic distances using dbscan. 
		"""
		self.pdbFileNames = [ x.replace(".pdb","") for x in pdbFileNames]
		self.pdbFileNames = [ x.replace(".kgs","") for x in self.pdbFileNames]
		# print('\n'.join(self.pdbFileNames))
		
		residueList = consistentResidueList(pdbFileNames)
		X = fullDataMatrix(pdbFileNames,residueList) 
		# printMatToFile(X,'shp2_distanceMatrix.txt') #Optional True to print full data matrix to file
		# X = backboneTorsionMatrix(pdbFileNames,residueList) 
		# printMatToFile(X,'shp2_torsionMatrix.txt') #Optional True to print full data matrix to file
		
		dbscan = cluster.DBSCAN(eps=1,min_samples=2).fit(X)
		labels = dbscan.labels_
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		print labels
		
class dataPCA:

	def __init__(self, pdbFileNames):
		"""
		Computes PCA on the pairwise-distance matrices across the pdb structures
		"""
		self.pdbFileNames = [ x.replace(".pdb","") for x in pdbFileNames]
		self.pdbFileNames = [ x.replace(".kgs","") for x in self.pdbFileNames]
		# print('\n'.join(self.pdbFileNames))
		residueList = consistentResidueList(pdbFileNames)
		X = fullDataMatrix(pdbFileNames,residueList)
		mu = np.mean(X, axis=0)
		n_components=3
		pca = PCA(n_components=n_components) #Select number of pcs, mle does automatic selection based on Thomas P. Minka: Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
		pca.fit(X)
		pca_score = pca.explained_variance_ratio_
		V = pca.components_
		print pca_score
		Xred = pca.transform(X,n_components)
		Xinv = pca.inverse_transform(Xred,n_components)
		# printMatToFile(Xred,"shp2_reducedData_3.txt")		
		# printMatToFile(Xinv,"shp2_recoveredData_3.txt")
		
		#Plot histograms to visualize density distribution
		# plt.hist(X,bins='auto')
		# plt.title("Distance distribution original data")
		# plt.savefig("shp2_originalDistribution.png")
		# plt.close()
		# 
		# plt.hist(Xinv,bins='auto')
		# plt.title("Distance distribution recovered data")
		# plt.savefig("shp2_recoveredDistribution.png")
		# plt.close()
		
		
		# x_pca_axis, y_pca_axis, z_pca_axis = V.T * pca_score / pca_score.min()
		# 
		# x_pca_axis, y_pca_axis, z_pca_axis = 3 * V.T
		# x_pca_plane = np.r_[x_pca_axis[:2], - x_pca_axis[1::-1]]
		# y_pca_plane = np.r_[y_pca_axis[:2], - y_pca_axis[1::-1]]
		# z_pca_plane = np.r_[z_pca_axis[:2], - z_pca_axis[1::-1]]
		# x_pca_plane.shape = (2, 2)
		# y_pca_plane.shape = (2, 2)
		# z_pca_plane.shape = (2, 2)
		# ax = createAxes(1,-40,-80)
		# ax.plot_surface(x_pca_plane, y_pca_plane, z_pca_plane)
		# ax.w_xaxis.set_ticklabels([])
		# ax.w_yaxis.set_ticklabels([])
		# ax.w_zaxis.set_ticklabels([])
		
	def createAxes(fig_num, elev, azim):
		fig = plt.figure(fig_num, figsize=(4, 3))
		plt.clf()
		ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=elev, azim=azim)
		return ax
		
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
	method='kmedoids'
	for i in range(len(sys.argv)):	
		if sys.argv[i]=='-kmedoids':
				method = 'kmedoids'
				del sys.argv[i]
				break;
		if sys.argv[i]=='-hierarchical':
				method = 'hierarchical'
				del sys.argv[i]
				break;
		if sys.argv[i]=='-dbscan':
				method = 'dbscan'
				del sys.argv[i]
				break;
		if sys.argv[i]=='-pca':
				method = 'pca'
				del sys.argv[i]
				break;
			
	print "Using "+method+" clustering"

	if(len(sys.argv)<=1):
		print "Usage:",sys.argv[0],"[-clusters <int>] [-kmedoids/-hierarchical] [-verbose] <list of pdb-files>"
		sys.exit(0)
	
	if method =='kmedoids':
		kc = KMedoidsClustering(sys.argv[1:],clusters)
		kc.printClusterInfo()
	elif method == 'hierarchical':
		hc = HierarchicalClustering(sys.argv[1:],clusters)
		hc.printClusterInfo()
	elif method == 'dbscan':
		db = DBscan(sys.argv[1:])
	elif method == 'pca':
		pca = dataPCA(sys.argv[1:])