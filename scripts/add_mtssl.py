#!/usr/bin/python
# coding: utf-8
"""
...
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

def scale(v, s):
    return [v[0]*s, v[1]*s, v[2]*s]

def add(v1, v2):
    return [ v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2] ]

def sub(v1, v2):
    return [ v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2] ]

def angle(v1, v2, v3):
    d21 = normalize( sub(v1,v2) )
    d23 = normalize( sub(v3,v2) )
    return math.acos( dot(d21, d23) )

def cross(v1,v2):
    return [ v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0] ]

def axisRotate(p_a, v_a, p, theta):
    """ Rotate p by the angle theta around the line defined by the point p_a and the vector v_a. Does not assume v_a is normalized. """
    #From http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    x,y,z = p
    a,b,c = p_a
    u,v,w = v_a
    t1 = 1-math.cos(theta) #Avoids repeated computations
    t2 = math.cos(theta)   #Avoids repeated computations
    t3 = math.sin(theta)   #Avoids repeated computations
    t4 = -u*x-v*y-w*z #Avoids repeated computations
    Rx = (a*(v*v+w*w)-u*(b*v+c*w+t4))*t1 + x*t2 + (-c*v+b*w-w*y+v*z)*t3
    Ry = (b*(u*u+w*w)-v*(a*u+c*w+t4))*t1 + y*t2 + ( c*u-a*w+w*x-u*z)*t3
    Rz = (c*(u*u+v*v)-w*(a*u+b*v+t4))*t1 + z*t2 + (-b*u+a*v-v*x+u*y)*t3
    return [Rx,Ry,Rz]



def putD(A, B, C, dist, theta, tau):
    #Initialize coordinate frame centered on C
    x = normalize( sub(B,C) )
    z = normalize( cross(sub(A,B),x) )
    y = cross(z,x)

    D = add(C, scale(x, dist)) #Place D at distance d from C along the C-B (-x) axis
    #print("ATOM   1000   SG CYS A  50    %8.3f%8.3f%8.3f  1.00  1.00          S"%(D[0],D[1],D[2]))
    D = axisRotate( C, z, D, -theta) #Rotate D around the z-axis
    #print("ATOM   1001   SG CYS A  50    %8.3f%8.3f%8.3f  1.00  1.00          S"%(D[0],D[1],D[2]))
    D = axisRotate( C, x, D, -tau)  #Rotate D around the x axis
    #print("ATOM   1002   SG CYS A  50    %8.3f%8.3f%8.3f  1.00  1.00          S"%(D[0],D[1],D[2]))
    return D


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



class Atom:
    """ Container class for PDB atoms """


    def __init__(self, id, name, alt, resn, chain, resi, pos, elem, hetatm):
        self.id = id
        self.name = name
        self.alt = alt
        self.resn = resn
        self.chain = chain
        self.resi = resi
        self.pos = pos
        self.elem = elem
        self.hetatm = hetatm

    def __str__(self):
        pre = "ATOM  "
        if self.hetatm: pre = "HETATM"
        return "%s%5d %4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00  1.00          %1s"%(pre,self.id, self.name,self.alt,self.resn,self.chain,self.resi,self.pos[0],self.pos[1],self.pos[2],self.elem)
        

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

#    def getSP(self):
#        if self.atomType in ["C3","N3+","N3","O3","S3+","S3","P3+"]:
#            return 3
#        if self.atomType in ["C2","Car","Cac","N2+","N2","Npl","Ng+","Ntr","O2","O3-","O2-","S2"]:
#            return 2
#        if self.atomType in ["C1","C1-","N1+","N1","O1+","O1"]:
#            return 1
#        raise RuntimeError("sp-hybridization not defined for atom type: "+self.atomType+" ("+str(self)+")")
#
#    def isDonor(self):
#        hasHNeighbor = any([n.elem=="H" for n in self.neighbors])
#        return hasHNeighbor and (self.elem=="O" or self.elem=="N")
#
#    def isAcceptor(self):
#        return self.elem=="O" or (self.elem=="N" and len(self.neighbors)<=2)


class PDBFile:
    """ A representation of a single-model PDB-file """

    def cleanAlternates(self):
        self.atoms = [a for a in self.atoms if a.alt in [" ","A"]]


    def getNearby(self, v, radius):
        irad = int(math.ceil(radius))
        ivx,ivy,ivz = int(v[0]+0.5), int(v[1]+0.5), int(v[2]+0.5)
        nearby = []

        for ix in range(ivx-irad, ivx+irad):
            for iy in range(ivy-irad, ivy+irad):
                for iz in range(ivz-irad, ivz+irad):
                    key = (ix,iy,iz)
                    if key in self.grid:
                        nearby+=[a for a in self.grid[key] if norm(sub(a.pos,v))<=radius]

        return nearby

#    def buildCovBonds(self):
#        """Use spatial hashing to update the `neighbors` field in each atom"""
#        for atom1 in self.atoms:
#            atom1.neighbors = []
#            #for atom2 in self.atoms: #All pairs
#                #if atom1==atom2: break #All pairs
#            for atom2 in self.getNearby(atom1.pos, atom1.covRadius()+2.0): #2.0 is the largest imaginable covalent radius for atom2
#                if atom1.id<=atom2.id: continue
#                dist = atom1.distance(atom2)
#                covsum = atom1.covRadius()+atom2.covRadius()
#                if dist<=covsum+0.4:
#                    atom1.neighbors.append(atom2)
#                    atom2.neighbors.append(atom1)

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

                a_id = int(line[6:11])
                a_name = line[11:16].strip()
                a_alt  = line[16]
                a_resn = line[17:20].strip()
                a_chain = line[21]
                a_resi = int(line[22:26])
                a_pos = [ float(line[30:38]), float(line[38:46]), float(line[46:54]) ]
                if len(line)>=78:
                    a_elem = line[76:78].strip()
                else:
                    a_elem = a_name[0]
                a_hetatm = line.startswith("HETATM")
                atom = Atom(a_id,a_name,a_alt,a_resn,a_chain,a_resi,a_pos,a_elem,a_hetatm)

                if a_alt in [" ","A"]:
                    self.atoms.append(atom)

            if line[0:6]=="ENDMDL" or line[0:5]=="MODEL":
                raise RuntimeError("Doesn't support multi-model pdbs")

        f.close()

        self.grid = defaultdict(list)
        for atom in self.atoms:
            self.grid[ (int(atom.pos[0]), int(atom.pos[1]), int(atom.pos[2])) ].append(atom)

        self.cleanAlternates()
        #self.buildCovBonds()


    def __repr__(self):
        return "PDBFile('"+self.file_name+"')"

    def __str__(self):
        return "\n".join([str(a) for a in self.atoms])

    def getAtom(self, resi, name):
        alist = [a for a in self.atoms if a.resi==resi and a.name==name]
        if len(alist)==1: return alist[0]
        return None

    def printLocalCBpos(self):
        resis = list(set([a.resi for a in self.atoms if a.resn!="GLY"]))
        avgs = [0,0,0]
        for resi in resis:
            N = self.getAtom(resi,"N")
            CA = self.getAtom(resi,"CA")
            C = self.getAtom(resi,"C")
            CB = self.getAtom(resi,"CB")

            #Place CB
            x = normalize( sub(C.pos, N.pos) )
            y = normalize( cross(sub(N.pos, CA.pos), x) )
            z = cross(x,y)
            CACB = sub(CB.pos, CA.pos)
            print("%.3f %.3f %.3f"%(dot(CACB,x),dot(CACB,y),dot(CACB,z)))
            
            avgs[0] += dot(CACB,x)
            avgs[1] += dot(CACB,y)
            avgs[2] += dot(CACB,z)
        avgs[0]/=len(resis)
        avgs[1]/=len(resis)
        avgs[2]/=len(resis)
        print "Average: "+str(avgs)

    def mutateToCYSstub(self, resi):
        """ Strip all side-chain atoms from residue and mutate to CYS """
        resiatoms = [a for a in self.atoms if a.resi==resi] # and not a.name in ["N", "CA", "C", "O", "H", "HA"] ]
        for a in resiatoms:
            if a.name in ["N", "CA", "C", "O", "H", "HA"]:
                a.resn="CYS"
            else:
                self.atoms.remove(a)
                self.grid[ (int(a.pos[0]), int(a.pos[1]), int(a.pos[2])) ].remove(a)

        #Place CB
        N = self.getAtom(resi,"N")
        CA = self.getAtom(resi,"CA")
        C = self.getAtom(resi,"C")
        x = normalize( sub(C.pos, N.pos) )
        y = normalize( cross(sub(N.pos, CA.pos), x) )
        z = cross(x,y)
        CBpos = add(CA.pos,add(scale(x,-0.036),add(scale(y,1.224),scale(z,-0.940))))
        CB = Atom(-1, "CB", " ", "CYS", CA.chain, CA.resi, CBpos, "C", False)
        self.atoms.append(CB)
        self.grid[ (int(CBpos[0]), int(CBpos[1]), int(CBpos[2])) ].append(CB)
            

    def addMTSSL(self, resi, rotamers):
        ret = []
        max_resi = max(self.atoms, lambda a: a.resi)[0].resi
        max_id   = max(self.atoms, lambda a: a.id)[0].id

        N = self.getAtom(resi,"N")
        CA = self.getAtom(resi,"CA")
        CB = self.getAtom(resi,"CB")
        vdw_factor = 0.9

        #Place remaining MTSSL probe
        SGpos = putD(N.pos,CA.pos,CB.pos, 1.76, 1.90, rotamers[0])
        ret.append(Atom(max_id+1, "SG", " ", "CYS", CA.chain, max_resi+1, SGpos, "S", True))
        #print ret[-1].name, [str(a) for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())], [str(a) for a in ret]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        S1pos = putD(CA.pos,CB.pos,SGpos, 2.00, 1.92, rotamers[1])
        ret.append(Atom(max_id+2, "S1", " ", "MTN", CA.chain, max_resi+1, S1pos, "S", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C4pos = putD(CB.pos,SGpos,S1pos, 1.76, 2.09, rotamers[2])
        ret.append(Atom(max_id+3, "C4", " ", "MTN", CA.chain, max_resi+1 , C4pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C3pos = putD(SGpos,S1pos,C4pos, 1.51, 1.88, rotamers[3])
        ret.append(Atom(max_id+4, "C3", " ", "MTN", CA.chain, max_resi+1 , C3pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C2pos = putD(S1pos,C4pos,C3pos, 1.34, 2.15, rotamers[4])
        ret.append(Atom(max_id+5, "C2", " ", "MTN", CA.chain, max_resi+1 , C2pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C1pos = putD(C4pos,C3pos,C2pos, 1.50, 1.92, 3.1415)
        ret.append(Atom(max_id+6, "C1", " ", "MTN", CA.chain, max_resi+1 , C1pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C8pos = putD(C3pos,C2pos,C1pos, 1.53, 1.88, -2.00)
        ret.append(Atom(max_id+7, "C8", " ", "MTN", CA.chain, max_resi+1 , C8pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C9pos = putD(C3pos,C2pos,C1pos, 1.53, 1.88, 2.00)
        ret.append(Atom(max_id+8, "C9", " ", "MTN", CA.chain, max_resi+1 , C9pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C5pos = putD(S1pos, C4pos, C3pos, 1.50, 2.16, rotamers[4]+3.1415)
        ret.append(Atom(max_id+9, "C5", " ", "MTN", CA.chain, max_resi+1 , C5pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C6pos = putD(C4pos, C3pos, C5pos, 1.52, 1.92, -1.07)
        ret.append(Atom(max_id+10, "C6", " ", "MTN", CA.chain, max_resi+1 , C6pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        C7pos = putD(C4pos, C3pos, C5pos, 1.52, 1.92,  1.07)
        ret.append(Atom(max_id+11, "C7", " ", "MTN", CA.chain, max_resi+1 , C7pos, "C", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        N1pos = putD(C4pos, C3pos, C5pos, 1.45, 1.76,  3.1415)
        ret.append(Atom(max_id+12, "N1", " ", "MTN", CA.chain, max_resi+1 , N1pos, "N", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        O1pos = putD(C3pos, C5pos, N1pos, 1.22, 2.14,  3.1415)
        ret.append(Atom(max_id+13, "O1", " ", "MTN", CA.chain, max_resi+1 , O1pos, "O", True))
        #print ret[-1].name, [str(a)+"\n" for a in self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius())]
        if self.getNearby(ret[-1].pos, vdw_factor*ret[-1].vdwRadius()): return None

        return ret


def pairPlusMinusDelta( R1, R2, d1, d2):
    ret = []
    for r1 in set([R1-d1, R1, R1+d1]):
        for r2 in set([R2-d2, R2, R2+d2]):
            ret.append( (r1,r2) )
    return ret
    


def genRotamersPMDelta(pdb):
    x1x2_combos = pairPlusMinusDelta( -76, -56, 10, 13)
    x1x2_combos+= pairPlusMinusDelta( -73, 173,  0,  0)
    x1x2_combos+= pairPlusMinusDelta(-170, -77,  3, 19)
    x1x2_combos+= pairPlusMinusDelta( 180, 180,  0,  0)
    x3_combos = [90, -90]
    x4x5_combos = [ (180,-77), (180,77), (-75,100), (-75,-8), (75,-100), (75,8) ]
    return [x1x2_combos, x3_combos, x4x5_combos]


if __name__ == "__main__":
    import sys

    pdb_file_name = sys.argv[1]
    resi = int(sys.argv[2])
    pdb = PDBFile(pdb_file_name)
    pdb.mutateToCYSstub(resi)
    
    mtssl_models = []
    rotamers = genRotamersPMDelta(pdb)
    for x1,x2 in rotamers[0]:
        for x3 in rotamers[1]:
            for x4,x5 in rotamers[2]:
                conf = pdb.addMTSSL(resi, [x1,x2,x3,x4,x5])
                if conf:
                    mtssl_models.append(conf)

    for i,conf in enumerate(mtssl_models):
        print "MODEL",i
        print str(pdb)
        print "\n".join([str(a) for a in conf])
        print "ENDMDL"
    



