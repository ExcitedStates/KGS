#!/usr/bin/python
# coding: utf-8
"""
Identify hybridization and hydrogen bonds in single-model PDB-files with hydrogens. 

A hydrogen bond is specified by a 5-tuple containing the acceptor base, acceptor, 
hydrogen, and donor atoms as well as h-bond energy. A donor atom is an oxygen or a 
nitrogen atom that has an attached hydrogen atom. An acceptor atom is any oxygen or 
a nitrogen with valence 2 or less. The acceptor base is the first neighbor of the 
acceptor.

Example:
    This script works on all types of molecules (protein, RNA/DNA, and ligands) 
    represented as PDB files. To print all h-bonds using this file as an external 
    module, construct a PDBFile object and call the getHydrogenBonds method:

        from hbfinder import PDBFile
        pdb = PDBFile("path/to/file.pdb")
        for aa,a,h,d,en in pdb.getHydrogenBonds():
            print("Bond %s - %s Energy: %.3f" % (str(a), str(d), en))


    If calling this script from the command-line type:

        ./hbfinder.py file.pdb

    This will print pymol commands to draw all hydrogen bonds.

References: 
    The hybridization algorithm is adapted from "Meng and Lewis. Determination of 
    Molecular Topology and Atomic Hybridization States from Heavy Atom Coordinates. 
    JCC. 1991". 

    The hydrogen bond energy is computed using the energy funtion in "Dahiyat, Gordon and Mayo.
    Automated design of the surface positions of protein helices. Protein Science, 1997."
"""

import math
import numpy as np
from collections import defaultdict

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
    d21 = normalize( sub(v1,v2) )
    d23 = normalize( sub(v3,v2) )
    return math.acos( dot(d21, d23) )

def cross(v1,v2):
    return [ v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] ]


#From TODO
covRadii = {"AC":1.88,"AG":1.59,"AL":1.35,"AM":1.51,"AS":1.21,"AU":1.50,"B":0.83 ,"BA":1.34,
            "BE":0.35,"BI":1.54,"BR":1.21,"C":0.68 ,"CA":0.99,"CD":1.69,"CE":1.83,"CL":0.99,
            "CO":1.33,"CR":1.35,"CS":1.67,"CU":1.52,"D":0.23 ,"DY":1.75,"ER":1.73,"EU":1.99,
            "F":0.64 ,"FE":1.34,"GA":1.22,"GD":1.79,"GE":1.17,"H":0.23 ,"HF":1.57,"HG":1.70,
            "HO":1.74,"I":1.40 ,"IN":1.63,"IR":1.32,"K":1.33 ,"LA":1.87,"LI":0.68,"LU":1.72,
            "MG":1.10,"MN":1.35,"MO":1.47,"N":0.68 ,"NA":0.97,"NB":1.48,"ND":1.81,"NI":1.50,
            "NP":1.55,"O":0.68 ,"OS":1.37,"P":1.05 ,"PA":1.61,"PB":1.54,"PD":1.50,"PM":1.80,
            "PO":1.68,"PR":1.82,"PT":1.50,"PU":1.53,"RA":1.90,"RB":1.47,"RE":1.35,"RH":1.45,
            "RU":1.40,"S":1.02 ,"SB":1.46,"SC":1.44,"SE":1.22,"SI":1.20,"SM":1.80,"SN":1.46,
            "SR":1.12,"TA":1.43,"TB":1.76,"TC":1.35,"TE":1.47,"TH":1.79,"TI":1.47,"TL":1.55,
            "TM":6.72,"U":1.58 ,"V":1.33 ,"W":1.37 ,"Y":1.78 ,"YB":1.94,"ZN":1.45,"ZR":1.56 }

cutoff=-1.0

class Atom:
    """ Container class for PDB atoms """

    def __init__(self, atom_string):
        self.id = int(atom_string[6:11])
        self.name = atom_string[11:16].strip()
        self.alt  = atom_string[16]
        self.resn = atom_string[17:20].strip()
        self.chain = atom_string[21]
        self.resi = int(atom_string[22:26])
        self.pos = [ float(atom_string[30:38]), float(atom_string[38:46]), float(atom_string[46:54]) ]
        if len(atom_string)>=78:
            self.elem = atom_string[76:78].strip()
        self.neighbors = []
        self.rings = 0
        self.atomType = ""
        self.hetatm = atom_string.startswith("HETATM")

    def __str__(self):
        return self.resn+str(self.resi)+"/"+self.name+"["+self.alt+","+self.elem+"]"

    def distance(self,a):
        return norm(sub(self.pos, a.pos))

    def vdwRadius(self):
        try:
            if self.elem=="C": return 1.7
            if self.elem=="N": return 1.625
            if self.elem=="O": return 1.480
            if self.elem=="S": return 1.782
            if self.elem=="H": return 1.0
        except AttributeError:
            pass
        return 1.8


    def covRadius(self):
        """The covalent radius of this atom."""
        try:
            return covRadii[self.elem.upper()]
        except AttributeError:
            pass
        raise RuntimeError("Unknown element for atom: "+str(self));

    def getSP(self):
        if self.atomType in ["C3","N3+","N3","O3","S3+","S3","P3+"]:
            return 3
        if self.atomType in ["C2","Car","Cac","N2+","N2","Npl","Ng+","Ntr","O2","O3-","O2-","S2"]:
            return 2
        if self.atomType in ["C1","C1-","N1+","N1","O1+","O1"]:
            return 1
        raise RuntimeError("sp-hybridization not defined for atom type: "+self.atomType+" ("+str(self)+")")

    def isDonor(self):
        hasHNeighbor = any([n.elem=="H" for n in self.neighbors])
        return hasHNeighbor and (self.elem=="O" or self.elem=="N")

    def isAcceptor(self):
        return self.elem=="O" or (self.elem=="N" and len(self.neighbors)<=2)


class PDBFile:
    """ A representation of a single-model PDB-file """

    def cleanAlternates(self):
        self.atoms = [a for a in self.atoms if a.alt in [" ","A"]]

    def cleanDehydroHOH(self):
        self.atoms = [a for a in self.atoms if not (a.resn=="HOH" and len(a.neighbors)==0)]


    def getNearby(self, v, radius):
        irad = int(math.ceil(radius))
        ivx,ivy,ivz = int(v[0]+0.5), int(v[1]+0.5), int(v[2]+0.5)
        nearby = []

        for ix in range(ivx-irad, ivx+irad):
            for iy in range(ivy-irad, ivy+irad):
                for iz in range(ivz-irad, ivz+irad):
                    key = (ix,iy,iz)
                    if key in self.grid:
                        nearby+=self.grid[key]

        return nearby

    def buildCovBonds(self):
        """Use spatial hashing to update the `neighbors` field in each atom"""
        for atom1 in self.atoms:
            atom1.neighbors = []
            #for atom2 in self.atoms: #All pairs
                #if atom1==atom2: break #All pairs
            for atom2 in self.getNearby(atom1.pos, atom1.covRadius()+2.0): #2.0 is the largest imaginable covalent radius for atom2
                if atom1.id<=atom2.id: continue
                dist = atom1.distance(atom2)
                covsum = atom1.covRadius()+atom2.covRadius()
                if dist<=covsum+0.4:
                    atom1.neighbors.append(atom2)
                    atom2.neighbors.append(atom1)

    def buildRingCounts(self):
        """Construct spanning tree to determine rings"""
        def ncaRing(T, v1, v2):
            v1path = []
            v = v1
            while v:
                v1path.append(v)
                v = T[v]
            v = v2
            v2path = []
            while not v in v1path:
                v2path.append(v)
                v = T[v]
            ring = v1path[0:v1path.index(v)+1] + v2path
            return ring
        T = {} #Associates an atom with its parent atom

        for root in self.atoms:
            if root in T: continue #Already explored
            T[root] = None
            fringe = [root]
            while fringe:
                a = fringe[0]
                del fringe[0]
                for n in a.neighbors:
                    if n in T and n == T[a]: continue # n is just parent of a
                    elif n in T and not (n in fringe): # There's a cycle
                        for r in ncaRing(T,a,n):
                            r.rings+=1
                    elif n not in fringe:
                        T[n] = a
                        fringe.append(n)

    def buildIDATM(self):
        """Follow the algorithm by Meng and Lewis 1991 to assign atom types. See 
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/idatm.html for a table of all types. """
        #Standard types for valences > 1
        for atom in self.atoms:
            if len(atom.neighbors)==4:
                if   atom.elem=="C": atom.atomType = "C3"
                elif atom.elem=="N": atom.atomType = "N3+"
                elif atom.elem=="P": 
                    freeOs = len([1 for n in atom.neighbors if n.elem=="O" and len(n.neighbors)==1])
                    if   freeOs>=2: atom.atomType = "Pac"
                    elif freeOs==1: atom.atomType = "Pox"
                    else:           atom.atomType = "P3+"
                elif atom.elem=="S":
                    freeOs = len([1 for n in atom.neighbors if n.elem=="O" and len(n.neighbors)==1])
                    if   freeOs>=3: atom.atomType = "Sac"
                    elif freeOs>=1: atom.atomType = "Son"
                    else:           atom.atomType = "S"

            elif len(atom.neighbors)==3:
                if   atom.elem=="C": 
                    freeOs = len([1 for n in atom.neighbors if n.elem=="O" and len(n.neighbors)==1])
                    atom.atomType = "C2" if freeOs<2 else "Cac"
                elif atom.elem=="N": 
                    freeOs = len([1 for n in atom.neighbors if n.elem=="O" and len(n.neighbors)==1])
                    atom.atomType = "Ntr" if freeOs>=2 else "Npl"
            elif len(atom.neighbors)==2:
                if   atom.elem=="C": atom.atomType = "C1"
                elif atom.elem=="N": atom.atomType = "N2"
                elif atom.elem=="O": atom.atomType = "O3"
                elif atom.elem=="S": atom.atomType = "S3"

        #Valence 1 atoms only
        for atom in self.atoms:
            if len(atom.neighbors)!=1: continue
            n = atom.neighbors[0]
            if   atom.elem=="C": 
                atom.atomType = "C1-" if n.elem=="O" and len(n.neighbors)==1 else "C1"
            elif atom.elem=="N": atom.atomType = "N1"
            elif atom.elem=="O": 
                if n.atomType=="Ntr" or n.atomType=="N1+" or n.atomType=="Cac":
                    atom.atomType = "O2-"
                elif n.atomType=="Pac" or n.atomType=="N3+" or n.atomType=="Sac" or n.atomType=="Son" or n.atomType=="Sxd":
                    atom.atomType="O3-"
                elif n.elem=="C": 
                    atom.atomType="O1+" if len(n.neighbors)==1 else "O2"
            elif atom.elem=="S": atom.atomType = "S2"
            elif atom.elem=="H": atom.atomType = "HC" if n.elem=="C" else "H"
            elif atom.elem=="D": atom.atomType = "DC" if n.elem=="C" else "D"


        #Check charge of nitrogens
        for atom in self.atoms:
            if atom.elem=="N" and atom.atomType!="N3+":
                if all([n.atomType in ["C3", "H", "D"] and n.rings==0 for n in atom.neighbors]):
                    atom.atomType="N3+"

            if atom.atomType=="C2":
                npls = [a for a in atom.neighbors if a.atomType in ["Npl", "Ng+"]]
                if len(npls)>=2 and all([n.rings==0 for n in npls]):
                    for n in npls:
                        n.atomType = "Ng+"




    def __init__(self, pdb):
        """Initialize a PDBFile object using a PDB-file or PDB-id. If pdb is 4
        characters long its assumed to be a PDB-id and the corresponding
        PDB-file will be downloaded and used. """

        self.file_name = pdb
        self.atoms = []

        if pdb.endswith(".gz"): 
            import gzip
            f = gzip.open(pdb, 'rb')
        else:
            f = open(pdb,'r')

        for line in f.readlines():
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                atom = Atom(line)
                if atom.alt in [" ","A"]:
                    self.atoms.append(atom)
            if line[0:6]=="ENDMDL" or line[0:5]=="MODEL":
                raise RuntimeError("Doesn't support multi-model pdbs")

        f.close()

        self.grid = defaultdict(list)
        for atom in self.atoms:
            self.grid[ (int(atom.pos[0]), int(atom.pos[1]), int(atom.pos[2])) ].append(atom)

        self.cleanAlternates()
        self.buildCovBonds()
        self.cleanDehydroHOH()
        self.buildRingCounts()
        self.buildIDATM()


    def hydrogenBondEnergy(self,aa, acceptor, hydrogen, donor):
        """Computes hydrogen bond energy following the terms described in
        "Dahiyat, Gordon and Mayo. Automated design of the surface positions of protein helices.
        Protein Science, 1997."
        If a hydrogen bond is not feasible, the value 1000 is returned.

        Args:
            aa (Atom): acceptor base
            acceptor (Atom): acceptor atom 
            hydrogen (Atom): hydrogen atom
            donor (Atom): donor atom

        Returns:
            float: 1000 is the geometric criteria for an hbond are not
            satisfied, or an h-bond energy otherwise
        """
        d0    = 8.0
        R0    = 2.8
        psi0  = 1.911 # 109.5degs
        R     = acceptor.distance(donor)
        theta = angle(donor.pos, hydrogen.pos, acceptor.pos)
        psi   = angle(hydrogen.pos, acceptor.pos, aa.pos)

        if R<2.6 or R>3.9: #R>3.2: #R>3.9
            return 1000
        if theta<(3.141592/2): 
            return 1000
        if donor.getSP()==3 and acceptor.getSP()==3 and psi-psi0>(3.141592/2):
            return 1000
        if donor.getSP()==3 and acceptor.getSP()==2 and psi<(3.141592/2):
            return 1000

        ratio       = R0/R
        ratio_pow2  = ratio*ratio
        ratio_pow4  = ratio_pow2*ratio_pow2
        ratio_pow8  = ratio_pow4*ratio_pow4
        ratio_pow10 = ratio_pow8*ratio_pow2
        ratio_pow12 = ratio_pow8*ratio_pow4

        energy_dist = d0 * (5*ratio_pow12 - 6*ratio_pow10)
        energy_angl = math.cos(theta)**2

        if(donor.getSP()==3 and acceptor.getSP()==3):
            return energy_dist * energy_angl * math.cos(psi-psi0)**2
        if(donor.getSP()==3 and acceptor.getSP()==2):
            return energy_dist * energy_angl * math.cos(psi)**2
        if(donor.getSP()==2 and acceptor.getSP()==3):
            return energy_dist * energy_angl * energy_angl
        if(donor.getSP()==2 and acceptor.getSP()==2):
            #compute phi
            n_d = normalize( cross(sub(donor.neighbors[0].pos, donor.pos), sub(donor.neighbors[1].pos,donor.pos)) )
            if( len(aa.neighbors) >=2 ):#Not necessarily the case
                n_a = normalize( cross(sub(aa.neighbors[0].pos, aa.pos), sub(aa.neighbors[1].pos,aa.pos)) )
            else:
                n_a = normalize( cross(sub(aa.pos,acceptor.pos), sub(hydrogen.pos,acceptor.pos)) )
            # n_a = normalize( cross(sub(aa.neighbors[0].pos, aa.pos), sub(aa.neighbors[1].pos,aa.pos)) )
            phi = math.acos( dot(n_d, n_a) )
            if phi>3.141592: phi = 3.141592-phi

            return energy_dist * energy_angl * math.cos(max(phi,psi))**2

        raise RuntimeError("Hybridization combination not implemented: "+str(donor.getSP())+" + "+str(acceptor.getSP()))

    
    def getHydrogenBonds(self, threshold=-1.0):
        """Goes through all relevant acceptor-donor pairs, finds base and hydrogen and collects
        it if the energy is below the threshold.

        Args:
            threshold (float): Energy threshold used to prune hydrogen bonds, default = -1.0

        Returns:
            list: Each element is a 5-tuple containing acceptor base (Atom), acceptor (Atom), 
            hydrogen (Atom), donor (Atom), and energy (float)
        """
        bonds = []
        for acceptor in self.atoms:
            if not acceptor.isAcceptor(): continue
            #for donor in self.atoms: # All-pairs
            for donor in self.getNearby(acceptor.pos, 3.9): #3.2
                if not donor.isDonor(): continue
                if acceptor==donor: continue
                if len(acceptor.neighbors)==0: continue
                aa = acceptor.neighbors[0]
                for hydrogen in donor.neighbors:
                    if not hydrogen.elem=="H": continue
                    if hydrogen.distance(acceptor) > 2.5: continue
                    try:
                        energy = self.hydrogenBondEnergy(aa, acceptor, hydrogen, donor)
                        if energy<threshold:
                            bonds.append( (aa,acceptor,hydrogen,donor, energy) )
                    except:
                        pass
        sorted(bonds, key = lambda x : min(x[1].id,x[2].id) )
        return bonds


    def burial(self, atom):
        """Computes the amount of atoms within a 8Å radius"""
        neighborhood = [a for a in self.getNearby(atom.pos, 8) if atom.distance(a)<8 ]
        return len(neighborhood)


    def __repr__(self):
        return "PDBFile('"+self.file_name+"')"



if __name__ == "__main__":
    import sys
    
    argv = sys.argv
    
    if "-energy" in argv:
        cutoff = sys.argv[sys.argv.index("-energy")+1]
        argv.remove("-energy")
        argv.remove(cutoff)
        print "Minimum h-bond energy "+cutoff
        cutoff = float(cutoff)
        
    for f in argv[1:]:
        print(f)
        pdb = PDBFile(f)

        for aa,a,h,d,energy in pdb.getHydrogenBonds(cutoff):
            # print(str(a.id)+" "+str(h.id)+" "+str(energy)) #KGS output
            # print("distance hbonds, ID "+str(a.id)+", ID "+str(h.id)) #Pymol output
            print(str(a.chain)+"/"+str(a.resi)+"/"+str(a.name)+" "+str(h.chain)+"/"+str(h.resi)+"/"+str(h.name)+" "+str(energy)) #KGS output
