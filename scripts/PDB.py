#!/usr/bin/python
"""
This file provides the classes PDB and Atom as a full PDB entry with all atoms
"""

class Atom:
    """ Container class for PDB atoms """

    def __init__(self, atom_string):
        self.id = int(atom_string[6:11])
        self.name = atom_string[11:16].strip()
        self.alt  = atom_string[16]
        self.resn = atom_string[17:20].strip()
        self.resn1 = oneLetterCode(self.resn)
        self.chain = atom_string[21]
        self.resi = int(atom_string[22:26])
        self.pos = [ float(atom_string[30:38]), float(atom_string[38:46]), float(atom_string[46:54]) ]
        self.occ = float(atom_string[54:60])
        self.tempFactor = float(atom_string[60:66])
        if (self.name[0] == "D"):#case of Deuterium
            self.name = "H"+self.name[1:]
        # if len(atom_string)>=78:
        #     self.elem = atom_string[76:78].strip()
        # else:
        self.elem = self.name[0]
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

    def getPDBline(self):
        line = "%-6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s \n" % ("HETATM" if self.hetatm else "ATOM", self.id, self.name, self.alt, self.resn, self.chain, self.resi, \
                self.pos[0], self.pos[1], self.pos[2], self.occ, self.tempFactor, self.elem)

        return line
    
class PDBFile:
    """ A representation of a single-model PDB-file """

    def cleanAlternates(self):
        self.atoms = [a for a in self.atoms if a.alt in [" ","A"]]

    def __init__(self, pdb):
        """Initialize a PDBFile object using a PDB-file or PDB-id."""

        self.file_name = pdb
        self.atoms = []

        f = open(pdb,'r')

        for line in f.readlines():
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                atom = Atom(line)
                if atom.alt in [" ","A"]:
                    self.atoms.append(atom)
            if line[0:6]=="ENDMDL" or line[0:5]=="MODEL":
                raise RuntimeError("Doesn't support multi-model pdbs")

        f.close()
        
        self.cleanAlternates()
        
    def getResi(self,atomID):
        for atom in self.atoms:
            if atom.id == atomID:
                return atom.resi;
            
    def getResn(self,atomID):
        for atom in self.atoms:
            if atom.id == atomID:
                return atom.resn;
            
    def getResnFromId(self,chain,resi):
        for atom in self.atoms:
            if atom.chain == chain and atom.resi == resi:
                return oneLetterCode(atom.resn);
            
    def getResn1(self,atomID):
        for atom in self.atoms:
            if atom.id == atomID:
                return atom.resn1;
            
    def getChain(self,atomID):
        for atom in self.atoms:
            if atom.id == atomID:
                return atom.chain;
            
    def getAtom(self,chain,resi,name):
        for atom in self.atoms:
            if atom.chain == chain and atom.resi == resi and atom.name==name:
                return atom;

    def writePDB(self,fname):
        f = open(fname,'w')

        for atom in self.atoms:
            f.write(atom.getPDBline())
            
        f.close()
        
def oneLetterCode(resn):
    return {
        'Arg': 'R',
        'ARG': 'R',
        'His': 'H',
        'HIS': 'H',
        'Hid': 'H',
        'HID': 'H',
        'Hie': 'H',
        'HIE': 'H',
        'Hip': 'H',
        'HIP': 'H',
        'Hse': 'H',
        'HSE': 'H',
        'Lys': 'K',
        'LYS': 'K',
        'Asp': 'D',
        'ASP': 'D',
        'Glu': 'E',
        'GLU': 'E',
        'Ser': 'S',
        'SER': 'S',
        'Thr': 'T',
        'THR': 'T',
        'Asn': 'N',
        'ASN': 'N',
        'Gln': 'Q',
        'GLN': 'Q',
        'Cys': 'C',
        'CYS': 'C',
        'Sec': 'U',
        'SEC': 'U',
        'Gly': 'G',
        'GLY': 'G',
        'Pro': 'P',
        'PRO': 'P',
        'Ala': 'A',
        'ALA': 'A',
        'Val': 'V',
        'VAL': 'V',
        'Ile': 'I',
        'ILE': 'I',
        'Leu': 'L',
        'LEU': 'L',
        'Met': 'M',
        'MET': 'M',
        'Phe': 'F',
        'PHE': 'F',
        'Tyr': 'Y',
        'TYR': 'Y',
        'Trp': 'W',
        'TRP': 'W',
    }.get(resn,'')