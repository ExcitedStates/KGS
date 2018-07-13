#!/usr/bin/python
# coding: utf-8
"""
Prepare input file for KGS from a PDB-file:
    1: Separates multi-model PDBs 
    2: Optionally fail if no hydrogens are present (and instruct to use reduce)
    3: Remove alt conformations
    4: Optionally remove waters
    5: Renumber atom ids
    6: Warn if steric collisions are present
    7: List covalent bonds for non-standard residues and HETATMs
    8: List disulphide bonds
    9: List hydrogen bonds

Example:
    This script works on all types of molecules (protein, RNA/DNA, and ligands) 
    represented as PDB files. 

        ./kgs_prepare.py 2N8B.pdb

    This will generate 20 files, 2N8B_{1..20}.kgs.pdb corresponding to each NMR 
    model and each will have 'REMARK KGS' records indicating hydrogen and disulphide 
    bonds. Had a ligand been present, records specifying rotatable bonds would also 
    be included.
"""

import math
import numpy as np
import random
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


covRadii = {"AC":1.88,"AG":1.59,"AL":1.35,"AM":1.51,"AS":1.21,"AU":1.50,"B" :0.83,"BA":1.34,
            "BE":0.35,"BI":1.54,"BR":1.21,"C" :0.68,"CA":0.99,"CD":1.69,"CE":1.83,"CL":0.99,
            "CO":1.33,"CR":1.35,"CS":1.67,"CU":1.52,"D" :0.23,"DY":1.75,"ER":1.73,"EU":1.99,
            "F" :0.64,"FE":1.34,"GA":1.22,"GD":1.79,"GE":1.17,"H" :0.23,"HF":1.57,"HG":1.70,
            "HO":1.74,"I" :1.40,"IN":1.63,"IR":1.32,"K" :1.33,"LA":1.87,"LI":0.68,"LU":1.72,
            "MG":1.10,"MN":1.35,"MO":1.47,"N" :0.68,"NA":0.97,"NB":1.48,"ND":1.81,"NI":1.50,
            "NP":1.55,"O" :0.68,"OS":1.37,"P" :1.05,"PA":1.61,"PB":1.54,"PD":1.50,"PM":1.80,
            "PO":1.68,"PR":1.82,"PT":1.50,"PU":1.53,"RA":1.90,"RB":1.47,"RE":1.35,"RH":1.45,
            "RU":1.40,"S" :1.02,"SB":1.46,"SC":1.44,"SE":1.22,"SI":1.20,"SM":1.80,"SN":1.46,
            "SR":1.12,"TA":1.43,"TB":1.76,"TC":1.35,"TE":1.47,"TH":1.79,"TI":1.47,"TL":1.55,
            "TM":6.72,"U" :1.58,"V" :1.33,"W" :1.37,"Y" :1.78,"YB":1.94,"ZN":1.45,"ZR":1.56 }

#From http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html#allatom, except for SE and UNKNOWN
vdwRadii = {
        "C" : 1.700, "N" : 1.625, "O" : 1.490, "S" : 1.782, "H" : 1.000, "P" : 1.871, 
        "F" : 1.560, "Cl": 1.735,"CL": 1.735,"Br": 1.978,"BR": 1.978, "I" : 2.094, "?" : 2.000 }

#Todo: implement ionic radii based on the coordination number
ionRadii = {
    "C" : 1 }

nulldistance={'C':1.908,'S':2.000}
welldepth={'C':-0.0860,'S':-0.2500}

verbose = False
writePml = False
writeIMOD = False
waters = True
ligands = True
hydrophobics = True
altloc= "A"
cutoff=-1.0
cutoffD=0.25 #Cutoff distance for hydrophobic interactions, sum of vdW + cutoffD
chainID="all"
model = -1

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
        self.occ = float(atom_string[54:60])
        self.tempFactor = float(atom_string[60:66])
        if (self.name[0] == "D"):#case of Deuterium
            self.name = "H"+self.name[1:]
        if len(atom_string.strip())>=78:
            self.elem = atom_string[76:78].strip()
            if self.isIon():
                self.elem=self.name
            else:
                self.elem = self.elem[0]
                # if self.elem == "D":
                #     self.elem = "H"
        else: #element column not provided in input
            name_nodigits = filter(lambda x: x.isalpha(), self.name)
            self.elem=name_nodigits
            if ~self.isIon():
                self.elem = name_nodigits[0]
       
        self.neighbors = []
        self.hBondNeighbors = []
        self.hydrophobicBondNeighbors=[]
        self.rings = 0
        self.atomType = ""
        self.hetatm = atom_string.startswith("HETATM")

    def __str__(self):
        return "//%s/%d/%s"%(self.chain, self.resi, self.name)

    def distance(self,a):
        return norm(sub(self.pos, a.pos))

    def vdwRadius(self):
        try:
            return vdwRadii[self.elem.upper()]
        except KeyError:
            return vdwRadii["?"]

    def covRadius(self):
        """The covalent radius of this atom."""
        try:
            return covRadii[self.elem.upper()]
        except AttributeError:
            pass
        except KeyError:
            RuntimeError("Unknown element for atom: "+str(self),self.name,self.atomType);
        raise RuntimeError("Unknown element for atom: "+str(self),self.name,self.elem);

    def welldepth(self):
        return welldepth[self.elem.upper()]

    def nulldistance(self):
        return nulldistance[self.elem.upper()]
    
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
        return self.elem=="O" or (self.elem=="N" and len(self.neighbors)<=2) # or self.elem == "S"

    def isHydrophobicAtom(self):
        #Limit hydrophobic interactions to C,S in non-polar residue side-chains
        return self.elem in ["C","S"] and self.resn in ["GLY","ALA","VAL","PRO","LEU","ILE","MET","TRP","PHE","CYS","TYR","GLN"] and not self.name=="C" and not self.name=="CA"
#        return self.elem in ["C","S"] and self.resn in ["GLY","ALA","VAL","PRO","LEU","ILE","MET","TRP","PHE","CYS","TYR","GLN"]:
    
    def getPDBline(self):
        if len(self.name)<=3: #Takes care of four-letter atom name alignment
            line = "%-6s%5d  %-3s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s \n" % ("HETATM" if self.hetatm else "ATOM", self.id, self.name, self.alt, self.resn, self.chain, self.resi, 
                self.pos[0], self.pos[1], self.pos[2], self.occ, self.tempFactor, self.elem)
        else:         
            line = "%-6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s \n" % ("HETATM" if self.hetatm else "ATOM", self.id, self.name, self.alt, self.resn, self.chain, self.resi, 
                    self.pos[0], self.pos[1], self.pos[2], self.occ, self.tempFactor, self.elem)

        return line
    
    def getCONECT(self):
        line=''
        if self.hetatm:
            line="CONECT%5d" % (self.id) # atom id, cols 7-11
            neighborCount=1
            for neighbor in self.neighbors:
                if not neighbor.hetatm: #Only covalent bonds for hetatms, other connections are revolute constraints in KGS
                    continue;
                if neighborCount > 4: #max. neighbors per CONECT entry (format requirement)
                    line += "\nCONECT%5d" % (self.id) #start new entry
                    neighborCount=1
                line += "%5d" % (neighbor.id) #neighbor ids, each 5-cols wide
                neighborCount += 1
            line += '\n' #new line
        return line
    
    def isIon(self):
        return self.elem in ["MG","ZN","K","FE"] #Todo: increase list for ions of interest
    
class PDBFile:
    """ A representation of a single-model PDB-file """

    def cleanDehydroHOH(self):
        dehydro_hoh = [a for a in self.atoms if      a.resn=="HOH" and len(a.neighbors)==0 ]
        self.atoms  = [a for a in self.atoms if not (a.resn=="HOH" and len(a.neighbors)==0)]
        if dehydro_hoh:
            print("%s: removed %d waters with no hydrogens"%(self.name, len(dehydro_hoh)))


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
            for atom2 in self.getNearby(atom1.pos, atom1.covRadius()+2.0): #2.0 is the largest imaginable covalent radius for atom2
                if atom1.id<=atom2.id: continue
                dist = atom1.distance(atom2)
                covsum = atom1.covRadius()+atom2.covRadius()
                if dist<=covsum+0.4:
                    atom1.neighbors.append(atom2)
                    atom2.neighbors.append(atom1)
        
    def rebuildCovBonds(self):
        """Use spatial hashing to update the `neighbors` field in each atom"""
        #Reser neighbor list
        for atom in self.atoms:
            atom.neighbors=[]
            
        for atom1 in self.atoms:
            for atom2 in self.getNearby(atom1.pos, atom1.covRadius()+2.0): #2.0 is the largest imaginable covalent radius for atom2
                if atom1.id<=atom2.id: continue
                dist = atom1.distance(atom2)
                covsum = atom1.covRadius()+atom2.covRadius()
                if dist<=covsum+0.4:
                    if (atom1.hetatm * atom2.hetatm): #Only covalent bonds between atom pairs or hetatm pairs
                        atom1.neighbors.append(atom2)
                        atom2.neighbors.append(atom1)
                        
    def buildHbondNeighbors(self):
        #Resetting neighbor information
        for atom in self.atoms:
            atom.hBondNeighbors=[]
            
        for aa,a,h,d,energy in self.getHydrogenBonds(cutoff):
            a.hBondNeighbors.append(h)
            h.hBondNeighbors.append(a)
            
    def buildHydrophobicBondNeighbors(self):
        for atom in self.atoms:
            atom.hydrophobicBondNeighbors=[]
        for c,s,energy in self.gethydrophobicBonds(cutoffD):
                c.hydrophobicBondNeighbors.append(s)  
                s.hydrophobicBondNeighbors.append(c)

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

    #Check metal ions (e.g. magnesium)
    def ionCoordination(self):
        for atom in self.atoms:
            if atom.isIon():
                for neighbor in atom.neighbors:
                    if neighbor.resn in ["A","G","C","U","T"]: #RNA/DNA
                        #......... TODO ................
                        # 1) identify neighbors for the ion
                        # 2) suitable distance constraints for ion interactions
                        # 3) Identify coordination for metal ions
                        pass
                    else:
                        print "MG ion close to protein:",atom.name,atom.id, neighbor.name,neighbor.id


    def __init__(self, name, atom_records):
        """ Initialize a PDBFile object using a PDB-file or PDB-id. """

        self.name = name
        self.atoms = atom_records
        self.constraints = []
        self.iModBonds = []
        
        self.checkHydrogens() #check if hydrogens are present
        self.checkAltConformations() #check and remove alt locs

        #Neighbor computations based on grid
        self.grid = defaultdict(list)
        for atom in self.atoms:
            self.grid[ (int(atom.pos[0]), int(atom.pos[1]), int(atom.pos[2])) ].append(atom)

        self.buildCovBonds() #Neighbors necessary to know which waters are present
        self.cleanDehydroHOH()  # check and potentially remove waters
 
        # #Redo grid after waters have been removed
        # self.grid = defaultdict(list)
        # for atom in self.atoms:
        #     self.grid[ (int(atom.pos[0]), int(atom.pos[1]), int(atom.pos[2])) ].append(atom)
            
        #Now we have all atoms in the model that we want --> no removal afterwards
        self.checkAtomIDs() #renumber atoms
        self.checkResidueSequence() #check residue sequence, identify odd numbering / gaps
        
        # self.rebuildCovBonds()
        self.buildRingCounts()
        self.buildIDATM()
        
        #Check ion neighborhood and coordination
        self.ionCoordination()


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

        if R<2.5 or R>3.9: #R>3.2: #R>3.9
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
        print "Identifying hbonds at energy cutoff ",threshold
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
                        # print "Angle: "+str(angle(donor.pos,hydrogen.pos,acceptor.pos)*180/math.pi)+", energy: "+str(energy)
                        if energy<threshold:
                            bonds.append( (aa,acceptor,hydrogen,donor, energy) )
                    except:
                        pass
        sorted(bonds, key = lambda x : min(x[1].id,x[2].id) )
        return bonds
    
    def gethydrophobicBonds(self,cutoffD=0.25):
        hybond=[]
        for c in self.atoms:
            if c.isHydrophobicAtom():
                for s in self.getNearby(c.pos,c.vdwRadius()+2.0+cutoffD):
                    if not s.isHydrophobicAtom(): continue;
                    if c.id>s.id: continue
                    if c.distance(s)<=(c.vdwRadius()+s.vdwRadius()+cutoffD) and not s.resi==c.resi:
                        energy=self.hydrophobicBondEnergy(c,s)
                        hybond.append((c,s,energy))
                        
        sorted(hybond, key = lambda x : min(x[0].id, x[1].id))
        return hybond
    
    def hydrophobicBondEnergy(self,c,s):
      
        '''Computes hydrophobic bond energy following the term described by Lennard Jones potential (as in KINARI)'''           
        if c in self.atoms and s in self.atoms: 
            if c.isHydrophobicAtom() and s.isHydrophobicAtom():
                        r=c.distance(s)
                        welldepth=math.sqrt(c.welldepth()*s.welldepth())
                        R=s.nulldistance()+c.nulldistance()
                        ratio=R/r
                        ratio2=ratio*ratio
                        ratio3=ratio2*ratio
                        ratio6=ratio3*ratio3
                        ratio12=ratio6*ratio6
                        return welldepth*(ratio12-2*ratio6)
        return 0

    def burial(self, atom):
        """Computes the amount of atoms within a 8Å radius"""
        neighborhood = [a for a in self.getNearby(atom.pos, 8) if atom.distance(a)<8 ]
        return len(neighborhood)


    def __repr__(self):
        return "PDBFile('"+self.name+"')"

    def checkHydrogens(self):
        hydrogens = [a for a in self.atoms if a.elem=="H"]
        if not hydrogens:
            print ("Warning: No hydrogens present in "+self.name+". Consider running `reduce` first.")
            # raise RuntimeError("No hydrogens present in "+self.name+". Consider running `reduce` first.")

    def checkAltConformations(self):
        if altloc=="all":
            return
        alt_confs =  [a for a in self.atoms if not (a.alt in [" ", altloc])]
        self.atoms = [a for a in self.atoms if      a.alt in [" ", altloc] ]
        for a in self.atoms:
            a.alt=" " #remove alternate conformation indicator
        if alt_confs:
            print(self.name+": Pruning "+str(len(alt_confs))+" atoms in alternative conformations")
        if alt_confs and verbose:
            for a in alt_confs:
                print("> removed alt conf "+str(a))

    def checkWaters(self):
        waters =     [a for a in self.atoms if a.resn=="HOH"]
        self.atoms = [a for a in self.atoms if a.resn!="HOH"]
        if waters:
            print(self.name+": Pruning "+str(len(waters))+" water atoms")
        if waters and verbose:
            for a in waters:
                print("> removed water "+str(a))


    def checkAtomIDs(self):
        for i,a in enumerate(self.atoms): a.id = i+1 #start at one

    def checkResidueSequence(self):
        currentResId = self.atoms[0].resi
        currentChain = self.atoms[0].chain
        resOffset = 0
        resList=[] #temporarily stores atoms of one residue
        inconsistencyFlag = 0
        for atom in self.atoms:
            if atom.hetatm: #Only need continuous sequence for residues
                continue
            if atom.chain != currentChain: #Reached a chain border, don't need to check for continous IDs
                currentChain = atom.chain #Set to new chain
                currentResId = atom.resi #Set to new first resi
                resList=[] #clear list
                resList.append((atom.name, atom.alt))
                if atom.alt != " ":  
                    resList.append((atom.name, " ")) #Fill list also with the "default" to correctly identify new residues
                if atom.alt == " ":
                    resList.append((atom.name, "A")) #Fill list also with the "default" to correctly identify new residues
                continue
            if atom.resi != currentResId or ((atom.name, atom.alt) in resList) : #New residue or wrong numbering
                # print atom.resi, atom.alt
                if atom.resi == currentResId + 1: #desired ID in contiguous sequence, go on
                    currentResId = atom.resi #move current residue one up
                    resList=[]
                    resList.append((atom.name, atom.alt))
                    if atom.alt != " ":  
                        resList.append((atom.name, " ")) #Fill list also with the "default" to correctly identify new residues
                    if atom.alt == " ":
                        resList.append((atom.name, "A")) #Fill list also with the "default" to correctly identify new residues
                    continue
                if atom.resi > currentResId + 1: #This is a gap
                    print "Chain "+str(currentChain)+": residue gap between "+str(currentResId)+" and "+str(atom.resi)
                    currentResId = atom.resi #move current residue across gap (user has to ensure correct sequence)
                    resList=[]
                    resList.append((atom.name, atom.alt))
                    if atom.alt != " ":  
                        resList.append((atom.name, " ")) #Fill list also with the "default" to correctly identify new residues
                    if atom.alt == " ":
                        resList.append((atom.name, "A")) #Fill list also with the "default" to correctly identify new residues
                    continue
                else: #weird numbering (e.g. 52 and 52A)
                    if ( (atom.name, atom.alt) in resList): #New residue (even if same wrong number, or same name)
                        if inconsistencyFlag == 0:
                            print "Inconsistent residue name at "+str(atom.resi)+atom.resn+" in chain "+str(currentChain)
                            inconsistencyFlag = 1
                        currentResId += 1
                        atom.resi = currentResId
                        resList=[]
                        resList.append((atom.name, atom.alt))
                        if atom.alt != " ":  
                            resList.append((atom.name, " ")) #Fill list also with the "default" to correctly identify new residues
                        if atom.alt == " ":
                            resList.append((atom.name, "A")) #Fill list also with the "default" to correctly identify new residues
                    else: #wrong numbering, in same residue
                        atom.resi = currentResId #adopt current ID
                        resList.append((atom.name, atom.alt))
                        if atom.alt != " ":  
                            resList.append((atom.name, " ")) #Fill list also with the "default" to correctly identify new residues
                        if atom.alt == " ":
                            resList.append((atom.name, "A")) #Fill list also with the "default" to correctly identify new residues
            else:      
                #Everything regular, just add name to list
                resList.append(atom.name) #keep track of atom names (each atom only once)
        
    def checkCollisions(self):
        serious_collisions = []
        collisions = []
        for atom1 in self.atoms:
            for atom2 in self.getNearby(atom1.pos, atom1.vdwRadius()+2.0): #2.0 is the largest vdw radius for atom2
                if atom1.id<=atom2.id: continue
                if abs(atom1.resi-atom2.resi)<=1 and atom1.chain==atom2.chain: continue
                dist = atom1.distance(atom2)
                if dist<0.7*(atom1.covRadius()+atom2.covRadius()):
                    serious_collisions.append( (atom1,atom2, dist) )
                    collisions.append( (atom1,atom2, dist) )
                elif dist<0.6*(atom1.vdwRadius()+atom2.vdwRadius()):
                    collisions.append( (atom1,atom2, dist) )
        if collisions or serious_collisions:
            print("%s: %d collisions, %d very serious"%(self.name, len(collisions), len(serious_collisions)))
            if verbose:
                for a1,a2,d in collisions:
                    print("> distance from %s to %s : %.2f"%(str(a1),str(a2),d))


    def checkCovalentBonds(self):
        odd_atoms = []
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="C" and not (len(a.neighbors) in [1,2,3,4])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="N" and not (len(a.neighbors) in [1,2,3,4])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="O" and not (len(a.neighbors) in [1,2])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="H" and not (len(a.neighbors) in [1])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="S" and not (len(a.neighbors) in [2,3,4,5,6])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem=="P" and not (len(a.neighbors) in [2,3,4,5])]
        if odd_atoms:
            print("%s: %d atoms with irregular number of covalent neighbors"%(self.name, len(odd_atoms)))
        if odd_atoms and verbose:
            for a in odd_atoms:
                print("> %s has %d covalent neighbors"%(str(a), len(a.neighbors)))


    def checkDisulphideBonds(self):
        sulphurs = [a for a in self.atoms if a.resn=="CYS" and a.elem=="S"]
        for s1 in sulphurs:
            s2 = [n for n in s1.neighbors if n.elem=="S"]
            if not s2:
                continue
            s2=s2[0]
            if s1.id>s2.id: continue
            self.constraints.append( "RevoluteConstraint %d %d %s %s"%(s1.id, s2.id, str(s1), str(s2)) )


    def checkHydrogenBonds(self):
        
        #Evaluate random h-bond set
        # maxID = max([atom.id for atom in self.atoms])
        # minID = min([atom.id for atom in self.atoms])
        # random.seed(418)
        # print maxID,minID
        # 
        # for aa,a,h,d,energy in self.getHydrogenBonds(cutoff):
        #     aID = random.randrange(minID,maxID)
        #     hID = random.randrange(minID,maxID)
        #     a = self.getAtom(aID)
        #     h = self.getAtom(hID)
        #     self.constraints.append( "RevoluteConstraint %d %d %s %s"%(a.id, h.id, str(a), str(h)) )
                
        # REAL VERSION
        for aa,a,h,d,energy in self.getHydrogenBonds(cutoff):
            self.constraints.append( "RevoluteConstraint %d %d %s %s"%(a.id, h.id, str(a), str(h)) )
            if writeIMOD:   #save h-bonds as stiff springs between donor and acceptor for usage in iMod
                firstNeighbor = d.neighbors[0]
                self.iModBonds.append( "%d %d %d"%(a.id,firstNeighbor.id, 10) )
                self.iModBonds.append( "%d %d %d"%(a.id, d.id, 10) )
                self.iModBonds.append( "%d %d %d"%(aa.id, d.id, 10) )

    def checkHydrophobicBonds(self):
        if hydrophobics:
            for a1,a2,energy in self.gethydrophobicBonds(cutoffD):
                self.constraints.append( "HydrophobicConstraint %d %d %s %s"%(a1.id, a2.id, str(a1), str(a2)) )

    def writePDB(self,fname):
        f = open(fname,'w')
        for con in self.constraints:
            f.write( "REMARK 555 %s\n"%(con) )

        for atom in self.atoms:
            f.write(atom.getPDBline())
        
        #Write constraints and CONECT records based on bond profile
        for atom in self.atoms:
            f.write(atom.getCONECT())
            
        f.close()
    
    def writePML(self,fname):
        f = open(fname,'w')
        weak=False;
        medium=False;
        strong=False;
        hydroph=False;
        for aa,a,h,d,energy in self.getHydrogenBonds(cutoff):
            if energy > -1.0:
                f.write( "distance hbonds_weak = id %s, id %s\n"%(a.id, h.id) )
                weak=True
            if energy <= -1.0 and energy > -3.0:
                f.write( "distance hbonds_medium = id %s, id %s\n"%(a.id, h.id) )
                medium=True
            if energy <= -3.0:
                f.write( "distance hbonds_strong = id %s, id %s\n"%(a.id, h.id) )
                strong=True
        for a1,a2,energy in self.gethydrophobicBonds(cutoffD):
            f.write("distance hydrophobics = id %s, id %s\n"%(a1.id,a2.id))
            hydroph=True
        if weak:
            f.write( "color white, hbonds_weak\n")
        if medium:
            f.write( "color yellow, hbonds_medium\n")
        if strong:
            f.write( "color red, hbonds_strong\n")
        if hydroph:
            f.write("color cyan, hydrophobics\n")
        f.write("hide labels\n")
        f.write("show cartoon\n")
        f.write("bg_color white\n")
        f.write("zoom\n")
        f.close()

    def getAtom(self, resi, name):
        ret = [a for a in self.atoms if a.resi==resi and a.name==name]
        if ret: return ret[0]
        return None

    def getAtom(self, atomId):
        ret = [a for a in self.atoms if a.id==atomId]
        if ret: return ret[0]
        return None
    
    def writeIMOD(self,fname):
        f = open(fname,'w')       
        for entry in self.iModBonds:
            f.write( "%s\n"%(entry) )
            
    def meanCoordination(self):
        self.buildHbondNeighbors()
        self.buildHydrophobicBondNeighbors()
        meanCoordination = 0.0
        for atom in self.atoms:
            meanCoordination += len(atom.neighbors)+len(atom.hBondNeighbors)+len(atom.hydrophobicBondNeighbors)
        return meanCoordination/len(self.atoms)

class MultiModel:
    def __init__(self, pdb_file,modelNum):
        """ Initialize a multi-model object using a PDB file name. """

        if pdb_file.endswith(".gz"): 
            import gzip
            f = gzip.open(pdb_file, 'rb')
        else:
            f = open(pdb_file,'r')

        models = []
        atoms = []
        remarks = []
        modelCount = 0;
        for line in f.readlines():
            if line[0:5]=="MODEL":
                if atoms: raise RuntimeError("A MODEL was not stopped (ENDMDL) before a new one started")
            elif line[0:6]=="ENDMDL":
                if not atoms: raise RuntimeError("Reached ENDMDL without reading any atoms")
                if (modelNum == -1): # -1 for all models, modelIndex for specific model only
                    models.append(atoms)
                if (modelCount == modelNum):
                    models.append(atoms)
                    break; #Only keeping this model
                modelCount += 1;
                atoms = []
            elif line[0:4] == "ATOM":
                atom = Atom(line)
                if chainID=="all" or (atom.chain==chainID): #only keep atom if in specified chains
                    atoms.append(atom)
            elif line[0:6] == "HETATM":
                atom = Atom(line)
                if chainID=="all" or (atom.chain==chainID): #only keep atom if in specified chains
                    if ligands and atom.resn != "HOH": #only add ligands if specified
                        atoms.append(atom)
                    if waters and atom.resn == "HOH": #only add waters if specified
                        atoms.append(atom)
            # elif line[0:6]=="REMARK":
            #     

        if atoms:
            models.append(atoms)

        f.close()
        base_name = pdb_file[0:pdb_file.rfind(".")] #finds last ".", either .pdb or .ent
        if len(models)==1:
            self.pdbs = [PDBFile(base_name, models[0])]
        else:
            self.pdbs = [PDBFile(base_name+"_"+str(i), m) for i,m in enumerate(models)]

    def rebuildCovalentBonds(self):
        for m in self.pdbs: m.rebuildCovBonds()
        
    def checkHydrogens(self):
        for m in self.pdbs: m.checkHydrogens()

    def checkAltConformations(self):
        for m in self.pdbs: m.checkAltConformations()
     
    def checkWaters(self):
        for m in self.pdbs: m.checkWaters()

    def checkAtomIDs(self):
        for m in self.pdbs: m.checkAtomIDs()

    def checkResidueSequences(self):
        for m in self.pdbs: m.checkResidueSequence()
        
    def checkCollisions(self):
        for m in self.pdbs: m.checkCollisions()

    def checkCovalentBonds(self):
        for m in self.pdbs: m.checkCovalentBonds()

    def checkDisulphideBonds(self):
        for m in self.pdbs: m.checkDisulphideBonds()

    def checkHydrogenBonds(self):
        for m in self.pdbs: m.checkHydrogenBonds()
        
    def checkHydrophobicBonds(self):
        for m in self.pdbs: m.checkHydrophobicBonds()

    def writePDBs(self, prefix):        
        for i,m in enumerate(self.pdbs):
            if len(self.pdbs)==1:
                m_prefix = "%s"%(prefix) if prefix else m.name
            else:
                m_prefix = "%s_%d"%(prefix,i) if prefix else m.name
            fname = m_prefix+".kgs.pdb"
            print("Writing "+m.name+" to "+fname)
            m.writePDB(fname)

    def writePMLs(self, prefix):
        for i,m in enumerate(self.pdbs): 
            if len(self.pdbs)==1:
                m_prefix = "%s"%(prefix) if prefix else m.name
            else:
                m_prefix = "%s_%d"%(prefix,i) if prefix else m.name
            fname = m_prefix+".kgs.pml"
            print("Writing constraints to "+fname)
            m.writePML(fname)
 
    def writeIMODs(self, prefix):
       for i,m in enumerate(self.pdbs): 
           if len(self.pdbs)==1:
               m_prefix = "%s"%(prefix) if prefix else m.name
           else:
               m_prefix = "%s_%d"%(prefix,i) if prefix else m.name
           fname = m_prefix+"_imod_ff.txt"
           print("Writing h-bond springs to "+fname)
           m.writeIMOD(fname)   

def printUsage(argv):
    print("Reads a PDB-file, prepares it for execution in KGS and stores the result")
    print("as a cleaned up PDB-file with REMARK records indicating constraints. The")
    print("output PDB will have the name '<prefix>_x.kgs.pdb' where the prefix will")
    print("be either the input-file name or whatever is specified by the -p option.")
    print("")
    print("usage: %s [options] <input>.pdb"%(argv[0]))
    print("where [options] can be either of the following")
    print("  -v            : verbose stdout")
    print("  -pre <string> : output-file prefix. Default is '<input>' ")
    print("  -pymol        : output a .pml file that displays constraints")
    print("  -noWaters     : remove all waters ")
    print("  -noLigands    : remove all ligands")
    print("  -noHydro      : no hydrophobic interactions")
    print("  -alt          : keeping specified alt loc (and empty alt loc) ")
    print("  -chain        : keeping only specified chain")
    print("  -energy       : energy cutoff for hydrogen bonds")
    print("  -cutoffD      : distance cutoff for hydrophobic interactions")
    print("  -imod         : save donor-acceptor connections for use as springs in iMod")
    print("  -model        : limit analysis to single model number in multi-model file")
    sys.exit(-1)

if __name__ == "__main__":
    import sys

    argv = sys.argv
    
    if "-h" in argv or "-help" in argv:
        printUsage(argv)

    if "-v" in argv:
        verbose = True
        argv.remove("-v")

    prefix = None
    if "-pre" in argv:
        i = argv.index("-pre")
        try:
            prefix = argv[i+1]
            if prefix[0]=="-": printUsage(argv)
            del argv[i+1]
        except IndexError:
            printUsage(argv)
        del argv[i]
        
    if "-pymol" in argv:
        writePml = True
        argv.remove("-pymol")
        
    if "-noLigands" in argv:
        ligands = False
        argv.remove("-noLigands")
        print "Removing all ligands"
        
    if "-noWaters" in argv:
        waters = False
        argv.remove("-noWaters")
        print "Removing all waters"

    if "-noHydro" in argv:
        hydrophobics = False
        argv.remove("-noHydro")
        print "No hydrophobics considered"
        
    if "-alt" in argv:
        altloc = sys.argv[sys.argv.index("-alt")+1]
        argv.remove("-alt")
        argv.remove(altloc)
        print "Keeping alternate location "+altloc

    if "-chain" in argv:
        chainID = sys.argv[sys.argv.index("-chain")+1]
        argv.remove("-chain")
        argv.remove(chainID)
        print "Keeping only chain "+chainID
        
    if "-energy" in argv:
        cutoff = sys.argv[sys.argv.index("-energy")+1]
        argv.remove("-energy")
        argv.remove(cutoff)
        print "Minimum h-bond energy "+cutoff
        cutoff = float(cutoff)
        
    if "-cutoffD" in argv:
        cutoffD = sys.argv[sys.argv.index("-cutoffD")+1]
        argv.remove("-cutoffD")
        argv.remove(cutoffD)
        print "Hydrophobic distance cutoff "+cutoffD
        cutoffD = float(cutoffD)
        
    if "-imod" in argv:
        writeIMOD = True
        argv.remove("-imod")
        
    if "-model" in argv:
        model = sys.argv[sys.argv.index("-model")+1]
        argv.remove("-model")
        argv.remove(model)
        print "Keeping only model",model
        model = int(model)

    #Only thing left of argv should be the program name and input file
    if not len(argv)==2:
        printUsage(argv)
    input_file = argv[1]

    models = MultiModel(input_file,model)
    models.checkCollisions()
    models.checkCovalentBonds()
    models.checkDisulphideBonds()
    models.checkHydrogenBonds()
    models.checkHydrophobicBonds()
    models.writePDBs(prefix)
    # models.pdbs[0].writePDB("init.kgs.pdb")
    if writePml:
        models.writePMLs(prefix)
    if writeIMOD:
        models.writeIMODs(prefix)
    

    # testVals= np.arange(-6,0,0.25)
    # for cutoff in testVals:
    #     print "Mean coordination of first model",models.pdbs[0].meanCoordination(),"cutoff",cutoff
