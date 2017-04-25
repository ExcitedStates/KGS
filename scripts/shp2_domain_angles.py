#! /usr/bin/python
# Copyright (c) 2017 Dominik Budday, dominik.budday@fau.de

'''

DESCRIPTION

Based on domain-separation in shp2, this script computes various angles between the domains,
and distances between defined residues to come up with coarse descriptives of the activation
and allosteric mechanism of SHP2.
The script is easily useable for other molecules when adapting the domains accordingly.
'''

import sys
import os
from pdb_structure import PDBFile
import math
from pdb_molecularmass import domainCOM
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']}) #Helvetica, Avant Gard, Computer Modern Sans serif
# from pylab import rcParams
# rcParams['figure.figsize'] = 6.5, 6.5*3/4

ftSize=14
mpl.rc('xtick', labelsize=ftSize) 
mpl.rc('ytick', labelsize=ftSize)

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def norm(v):
    return math.sqrt(dot(v,v))

def normalize(v):
    l = norm(v)
    return [v[0]/l, v[1]/l, v[2]/l]

def sub(v1, v2):
    return [ v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2] ]

def angle(v1, v2, v3):
    #Angle between vectors (2->1, 2->3)
    d21 = normalize( sub(v1,v2) )
    d23 = normalize( sub(v3,v2) )
    return math.acos( dot(d21, d23) )

def cross(v1,v2):
    return [ v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] ]

def shp2DomainAngles(pdbStructure,availableResis):
    
    #3 Domains in SHP2 --> Adapt the residue list for useage with other proteins
    ptpRange = range(220,530)
    csh2Range = range(112,219)
    nsh2Range = range(1,111)
    
    #Truncate lists to residues available in all structures
    ptpRange = [x for x in ptpRange if x in availableResis]
    csh2Range = [x for x in csh2Range if x in availableResis]
    nsh2Range = [x for x in nsh2Range if x in availableResis]
    
    #Compute the Center of Mass for each domain
    ptpAtoms = []
    for resi in ptpRange:
        ptpAtoms.extend(pdbStructure.getAtomsInResi(resi))
    x,y,z = domainCOM(ptpAtoms)
    ptpCOM = [x,y,z]
    # print "PTP COM: "+str(x)+", "+str(y)+", "+str(z)

    csh2Atoms = []
    for resi in csh2Range:
        csh2Atoms.extend(pdbStructure.getAtomsInResi(resi))
    x,y,z = domainCOM(csh2Atoms)
    csh2COM = [x,y,z]
    # print "CSH2 COM: "+str(x)+", "+str(y)+", "+str(z)
    
    nsh2Atoms = []
    for resi in nsh2Range:
        nsh2Atoms.extend(pdbStructure.getAtomsInResi(resi))
    x,y,z = domainCOM(nsh2Atoms)
    nsh2COM = [x,y,z]
    # print "NSH2 COM: "+str(x)+", "+str(y)+", "+str(z)

    #Full COM
    allAtoms = nsh2Atoms
    allAtoms.extend(csh2Atoms)
    allAtoms.extend(ptpAtoms)
    x,y,z = domainCOM(allAtoms)
    fullCOM = [x,y,z]
    # print "Overall COM: "+str(x)+", "+str(y)+", "+str(z)
    
    #Compute the angles between the centers
    alpha = angle(ptpCOM, fullCOM, nsh2COM)
    beta = angle(ptpCOM, fullCOM, csh2COM)
    delta = angle(nsh2COM, fullCOM, csh2COM)
    return alpha, beta, delta

def shp2Distances(pdbStructure):
    pass
    # To be done
    # distancePairs=[[111,253],[114,249],[111,216],[111,217]]

def consistentResidueList(pdbFileNames):
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

def shp2PrincipalAxesAngles(pdbStructure,availableResis):
    
    #3 Domains in SHP2 --> Adapt the residue list for useage with other proteins
    ptpRange = range(220,530)
    csh2Range = range(112,219)
    nsh2Range = range(1,111)
    
    #Truncate lists to residues available in all structures
    ptpRange = [x for x in ptpRange if x in availableResis]
    csh2Range = [x for x in csh2Range if x in availableResis]
    nsh2Range = [x for x in nsh2Range if x in availableResis]
    
    #Compute the Center of Mass for each domain
    ptpAtoms = []
    for resi in ptpRange:
        ptpAtoms.extend(pdbStructure.getAtomsInResi(resi))
    ptpCoordMat = [[a.x, a.y, a.z] for a in ptpAtoms]

    pca = PCA(n_components=1) #Select number of pcs, mle does automatic selection based on Thomas P. Minka: Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    pca.fit(ptpCoordMat)
    pca_score = pca.explained_variance_ratio_
    ptpV = pca.components_
    print ptpV
    x,y,z = domainCOM(ptpAtoms)
    ptpCOM = [x,y,z]
    print "PTP COM: "+str(x)+", "+str(y)+", "+str(z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ptpCoordMat[:][0],ptpCoordMat[:][1],ptpCoordMat[:][2])
    
    csh2Atoms = []
    for resi in csh2Range:
        csh2Atoms.extend(pdbStructure.getAtomsInResi(resi))
    csh2CoordMat = [[a.x, a.y, a.z] for a in csh2Atoms]
    print csh2CoordMat[0][:]
    pca = PCA(n_components=1) #Select number of pcs, mle does automatic selection based on Thomas P. Minka: Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    pca.fit(csh2CoordMat)
    pca_score = pca.explained_variance_ratio_
    V = pca.components_
    print V
    x,y,z = domainCOM(csh2Atoms)
    csh2COM = [x,y,z]
    print "CSH2 COM: "+str(x)+", "+str(y)+", "+str(z)

    nsh2Atoms = []
    for resi in nsh2Range:
        nsh2Atoms.extend(pdbStructure.getAtomsInResi(resi))
    nsh2CoordMat = [[a.x, a.y, a.z] for a in nsh2Atoms]
    print nsh2CoordMat[0][:]
    pca = PCA(n_components=1) #Select number of pcs, mle does automatic selection based on Thomas P. Minka: Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
    pca.fit(nsh2CoordMat)
    pca_score = pca.explained_variance_ratio_
    V = pca.components_
    print V
    x,y,z = domainCOM(nsh2Atoms)
    nsh2COM = [x,y,z]
    print "NSH2 COM: "+str(x)+", "+str(y)+", "+str(z)
    
    return x,y

if __name__ == "__main__":
    if len(sys.argv)<2:
        print "Usage: "+sys.argv[0]+" <pdb-files in a row>"
        sys.exit(1)
    
    deltaSH2=[]
    ptp_NSH2=[]
    ptp_CSH2=[]
    names=[]
    indexes = range(0,len(sys.argv)-1)
    pdbFileNames = sys.argv[1:]
    availableResidues = consistentResidueList(pdbFileNames)
    
    for pdb_file in pdbFileNames:
        struc = PDBFile(pdb_file) 
    
        name = os.path.basename(pdb_file).replace(".pdb","")
        name = name.replace(".kgs","")
        
        alpha, beta, delta = shp2DomainAngles(struc,availableResidues)
        # print alpha*180/math.pi, beta*180/math.pi, delta*180/math.pi
        
        names.append(name)
        ptp_NSH2.append(alpha*180/math.pi)
        ptp_CSH2.append(beta*180/math.pi)
        deltaSH2.append(delta*180/math.pi)
        
        #Determine principal axes and angles between them
        # phi, psi = shp2PrincipalAxesAngles(struc,availableResidues)
        # sys.exit(1)
     
    #Sort them by some thing: in this case delta angles   
    sortedIndices = [i[0] for i in sorted(enumerate(deltaSH2), key=lambda x:x[1])]
    # print sorted(enumerate(deltaSH2), key=lambda x:x[1])
    ptp_NSH2 = [ptp_NSH2[x] for x in sortedIndices]
    ptp_CSH2 = [ptp_CSH2[x] for x in sortedIndices]
    deltaSH2 = [deltaSH2[x] for x in sortedIndices]
    names = [names[x] for x in sortedIndices]
    
    fig = plt.figure()
    ax = fig.add_axes([0.2,0.1,0.6,0.85])
    # ax.plot(names,{ptp_NSH2,deltaSH2})
    plt.plot(ptp_NSH2,indexes,'ro',ptp_CSH2,indexes,'bo',deltaSH2,indexes,'go')
    ax.set_ylim(-0.5,len(sys.argv)-1.5)
    
    indexesLeft = indexes[0::2]
    namesLeft = names[0::2]
    indexesRight = indexes[1::2]
    namesRight = names[1::2]

    ax.set_yticks(indexesLeft)
    ax.set_yticklabels(namesLeft,fontsize=ftSize)
    
    ax2 = ax.twinx()
    ax2.set_ylim(-0.5,len(sys.argv)-1.5)
    ax2.set_yticks(indexesRight)
    ax2.set_yticklabels(namesRight,fontsize=ftSize)
        
    ax.set_xlabel('Domain angle [$^\circ$]',fontsize=ftSize)
    ax.grid(True)
    ax2.grid(True)
    # plt.xticks(indexes,names,rotation='vertical')
    # plt.margins(0.2)
    # plt.subplots_adjust(bottom=0.2)
    plt.savefig("shp2_domain_angles.png", transparent='True')
    
    orig_stdout = sys.stdout
    fout = file('shp2_domain_angles.txt','w')
    sys.stdout = fout
    print "Name PTP_NSH2_angle PTP_CSH2_angle SH2_delta_angle"
    for i in range(len(names)):
        print names[i]+" "+str(ptp_NSH2[i])+" "+str(ptp_CSH2[i])+" "+str(deltaSH2[i])
       
    fout.close()
    sys.stdout = orig_stdout