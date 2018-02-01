#!/usr/bin/env python2
# coding: utf-8
"""
Classes for representing and manipulating pdb-structures.

Example:
    This library works on all types of molecules (protein, RNA/DNA, and ligands)
    represented as PDB files. To print all h-bonds using this file as an external 
    module, construct a PDBFile object and call the get_hydrogen_bonds method:

        from pdb_structure import *
        pdb = pdb_structure.PDBFile("path/to/file.pdb")
        for aa, a, h, d, en in pdb.compute_hydrogen_bonds():
            print("Bond %s - %s Energy: %.3f" % (str(a), str(d), en))

References:
    The hybridization algorithm is adapted from "Meng and Lewis. Determination of 
    Molecular Topology and Atomic Hybridization States from Heavy Atom Coordinates. 
    JCC. 1991". 

    The hydrogen bond energy is computed using the energy funtion in "Dahiyat, Gordon and Mayo.
    Automated design of the surface positions of protein helices. Protein Science, 1997."
"""
import math
from collections import defaultdict

verbose = False


def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def norm(v):
    return math.sqrt(dot(v, v))


def normalize(v):
    l = norm(v)
    return [v[0]/l, v[1]/l, v[2]/l]


def sub(v1, v2):
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]


def angle(v1, v2, v3):
    d21 = normalize(sub(v1, v2))
    d23 = normalize(sub(v3, v2))
    return math.acos(dot(d21, d23))


def cross(v1, v2):
    return [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]


cov_radii = {"AC": 1.88, "AG": 1.59, "AL": 1.35, "AM": 1.51, "AS": 1.21, "AU": 1.50, "B":  0.83, "BA": 1.34,
             "BE": 0.35, "BI": 1.54, "BR": 1.21, "C":  0.68, "CA": 0.99, "CD": 1.69, "CE": 1.83, "CL": 0.99,
             "CO": 1.33, "CR": 1.35, "CS": 1.67, "CU": 1.52, "D":  0.23, "DY": 1.75, "ER": 1.73, "EU": 1.99,
             "F":  0.64, "FE": 1.34, "GA": 1.22, "GD": 1.79, "GE": 1.17, "H":  0.23, "HF": 1.57, "HG": 1.70,
             "HO": 1.74, "I":  1.40, "IN": 1.63, "IR": 1.32, "K":  1.33, "LA": 1.87, "LI": 0.68, "LU": 1.72,
             "MG": 1.10, "MN": 1.35, "MO": 1.47, "N":  0.68, "NA": 0.97, "NB": 1.48, "ND": 1.81, "NI": 1.50,
             "NP": 1.55, "O":  0.68, "OS": 1.37, "P":  1.05, "PA": 1.61, "PB": 1.54, "PD": 1.50, "PM": 1.80,
             "PO": 1.68, "PR": 1.82, "PT": 1.50, "PU": 1.53, "RA": 1.90, "RB": 1.47, "RE": 1.35, "RH": 1.45,
             "RU": 1.40, "S":  1.02, "SB": 1.46, "SC": 1.44, "SE": 1.22, "SI": 1.20, "SM": 1.80, "SN": 1.46,
             "SR": 1.12, "TA": 1.43, "TB": 1.76, "TC": 1.35, "TE": 1.47, "TH": 1.79, "TI": 1.47, "TL": 1.55,
             "TM": 6.72, "U":  1.58, "V":  1.33, "W":  1.37, "Y":  1.78, "YB": 1.94, "ZN": 1.45, "ZR": 1.56
             }


# From http://www.cgl.ucsf.edu/chimera/docs/Users_guide/midas/vdwtables.html# allatom, except for SE and UNKNOWN
vdw_radii = {"C": 1.700, "N":  1.625, "O":  1.490, "S":  1.782, "H":  1.000, "P": 1.871,
             "F": 1.560, "Cl": 1.735, "CL": 1.735, "Br": 1.978, "BR": 1.978, "I": 2.094, "?": 2.000
             }


class Atom:
    """ Container class for PDB atoms """

    def __init__(self, atom_string):
        self.id = int(atom_string[6:11])
        self.name = atom_string[11:16].strip()
        self.alt = atom_string[16]
        self.resn = atom_string[17:20].strip()
        self.chain = atom_string[21]
        self.resi = int(atom_string[22:26])
        self.pos = [float(atom_string[30:38]), float(atom_string[38:46]), float(atom_string[46:54])]
        self.occ = float(atom_string[54:60])
        self.temp_factor = float(atom_string[60:66])
        if self.name[0] == "D":  # case of Deuterium
            self.name = "H"+self.name[1:]
        if len(atom_string.strip()) >= 78:
            self.elem = atom_string[76:78].strip()
            # self.elem = self.elem[0]
            if self.elem == "D":
                self.elem = "H"
        else:
            name_nodigits = filter(lambda x: x.isalpha(), self.name)
            self.elem = name_nodigits[0]
        self.neighbors = []
        self.rings = 0
        self.atom_type = ""
        self.hetatm = atom_string.startswith("HETATM")

    def __str__(self):
        return "//%s/%d/%s" % (self.chain, self.resi, self.name)

    def __getitem__(self, key):
        return self.pos[key]

    def distance(self, a):
        return norm(sub(self, a))

    def vdwradius(self):
        try:
            return vdw_radii[self.elem.upper()]
        except KeyError:
            return vdw_radii["?"]

    def covradius(self):
        """The covalent radius of this atom."""
        try:
            return cov_radii[self.elem.upper()]
        except AttributeError:
            pass
        except KeyError:
            RuntimeError("Unknown element ("+self.elem+") for atom: "+str(self))
        raise RuntimeError("Unknown element ("+self.elem+") for atom: "+str(self))

    def get_sp(self):
        if self.atom_type in ["C3", "N3+", "N3", "O3", "S3+", "S3", "P3+"]:
            return 3
        if self.atom_type in ["C2", "Car", "Cac", "N2+", "N2", "Npl", "Ng+", "Ntr", "O2", "O3-", "O2-", "S2"]:
            return 2
        if self.atom_type in ["C1", "C1-", "N1+", "N1", "O1+", "O1"]:
            return 1
        raise RuntimeError("sp-hybridization not defined for atom type: "+self.atom_type+" ("+str(self)+")")

    def is_donor(self):
        has_h_neighbor = any([n.elem == "H" for n in self.neighbors])
        return has_h_neighbor and (self.elem == "O" or self.elem == "N")

    def is_acceptor(self):
        return self.elem == "O" or (self.elem == "N" and len(self.neighbors) <= 2)

    def pdb_str(self):
        if len(self.name) <= 3:  # Takes care of four-letter atom name alignment
            line = "%-6s%5d  %-3s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s \n" % \
                   ("HETATM" if self.hetatm else "ATOM", self.id, self.name, self.alt, self.resn, self.chain, self.resi,
                    self.pos[0], self.pos[1], self.pos[2], self.occ, self.temp_factor, self.elem)
        else:         
            line = "%-6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s \n" % \
                   ("HETATM" if self.hetatm else "ATOM", self.id, self.name, self.alt, self.resn, self.chain, self.resi,
                    self.pos[0], self.pos[1], self.pos[2], self.occ, self.temp_factor, self.elem)

        return line
    
    def get_conect(self):
        line = ''
        if self.hetatm:
            line = "CONECT%5d" % self.id  # atom id, cols 7-11
            neighbor_count = 1
            for neighbor in self.neighbors:
                if not neighbor.hetatm:  # Only covalent bonds for hetatms
                    continue
                if neighbor_count > 4:  # max. neighbors per CONECT entry (format requirement)
                    line += "\nCONECT%5d" % self.id  # start new entry
                    neighbor_count = 1
                line += "%5d" % neighbor. id  # neighbor ids, each 5-cols wide
                neighbor_count += 1
            line += '\n'  # new line
        return line


class PDBModel:
    """ A representation of a single-model PDB-file """

    def clean_dehydro_hoh(self):
        dehydro_hoh = [a for a in self.atoms if a.resn == "HOH" and len(a.neighbors) == 0]
        self.atoms = [a for a in self.atoms if not (a.resn == "HOH" and len(a.neighbors) == 0)]
        if dehydro_hoh:
            print("%s: removed %d waters with no hydrogens" % (self.name, len(dehydro_hoh)))

    def get_nearby(self, v, radius):
        """
        Locate all atoms within `radius` of `v`.
        :param v: A list of coordinates, or an `Atom`
        :param radius: A number
        :return: A list of atoms near `v`
        """
        return [a for a in self.get_approx_nearby(v, radius) if a.distance(v) <= radius]

    def get_approx_nearby(self, v, radius):
        """
        Locate all atoms within `radius` of `v`.
        :param v: A list of coordinates, or an `Atom`
        :param radius: A number
        :return: A list of atoms near `v`
        """
        irad = int(math.ceil(radius))
        ivx, ivy, ivz = int(v[0]+0.5), int(v[1]+0.5), int(v[2]+0.5)
        nearby = []

        for ix in range(ivx-irad, ivx+irad):
            for iy in range(ivy-irad, ivy+irad):
                for iz in range(ivz-irad, ivz+irad):
                    key = (ix, iy, iz)
                    if key in self.grid:
                        nearby += self.grid[key]

        return [a for a in nearby if a.distance(v) <= radius]

    def build_cov_bonds(self):
        """Use spatial hashing to update the `neighbors` field in each atom"""
        for atom1 in self.atoms:
            # 2.0 is the largest imaginable covalent radius for atom2
            for atom2 in self.get_approx_nearby(atom1.pos, atom1.covradius()+2.0):
                if atom1.id <= atom2.id:
                    continue
                dist = atom1.distance(atom2)
                covsum = atom1.covradius()+atom2.covradius()
                if dist <= covsum+0.4:
                    atom1.neighbors.append(atom2)
                    atom2.neighbors.append(atom1)
        
    def rebuild_cov_bonds(self):
        """Use spatial hashing to update the `neighbors` field in each atom"""
        # Reser neighbor list
        for atom in self.atoms:
            atom.neighbors = []
            
        for atom1 in self.atoms:
            # 2.0 is the largest imaginable covalent radius for atom2
            for atom2 in self.get_approx_nearby(atom1.pos, atom1.covradius()+2.0):
                if atom1.id <= atom2.id:
                    continue
                dist = atom1.distance(atom2)
                covsum = atom1.covradius()+atom2.covradius()

                if dist <= covsum+0.4:
                    if atom1.hetatm * atom2.hetatm:  # Only covalent bonds between atom pairs or hetatm pairs
                        atom1.neighbors.append(atom2)
                        atom2.neighbors.append(atom1)         

    def build_ring_counts(self):
        """Construct spanning tree to determine rings"""
        def nca_ring(t, v1, v2):
            v1path = []
            v = v1
            while v:
                v1path.append(v)
                v = t[v]
            v = v2
            v2path = []
            while v not in v1path:
                v2path.append(v)
                v = t[v]
            ring = v1path[0:v1path.index(v)+1] + v2path
            return ring

        parent_dict = {}  # Associates an atom with its parent atom

        for root in self.atoms:
            if root in parent_dict:
                continue  # Already explored
            parent_dict[root] = None
            fringe = [root]
            while fringe:
                a = fringe[0]
                del fringe[0]
                for n in a.neighbors:
                    if n in parent_dict and n == parent_dict[a]:
                        continue  # n is just parent of a
                    elif n in parent_dict and not (n in fringe):  # There's a cycle
                        for r in nca_ring(parent_dict, a, n):
                            r.rings += 1
                    elif n not in fringe:
                        parent_dict[n] = a
                        fringe.append(n)

    def build_idatm(self):
        """Follow the algorithm by Meng and Lewis 1991 to assign atom types. See 
        https://www.cgl.ucsf.edu/chimera/docs/Users_guide/idatm.html for a table of all types. """
        # Standard types for valences > 1
        for atom in self.atoms:
            if len(atom.neighbors) == 4:
                if atom.elem == "C":
                    atom.atom_type = "C3"
                elif atom.elem == "N":
                    atom.atom_type = "N3+"
                elif atom.elem == "P":
                    free_os = len([1 for n in atom.neighbors if n.elem == "O" and len(n.neighbors) == 1])
                    if free_os >= 2:
                        atom.atom_type = "Pac"
                    elif free_os == 1:
                        atom.atom_type = "Pox"
                    else:
                        atom.atom_type = "P3+"
                elif atom.elem == "S":
                    free_os = len([1 for n in atom.neighbors if n.elem == "O" and len(n.neighbors) == 1])
                    if free_os >= 3:
                        atom.atom_type = "Sac"
                    elif free_os >= 1:
                        atom.atom_type = "Son"
                    else:
                        atom.atom_type = "S"

            elif len(atom.neighbors) == 3:
                if atom.elem == "C":
                    free_os = len([1 for n in atom.neighbors if n.elem == "O" and len(n.neighbors) == 1])
                    atom.atom_type = "C2" if free_os < 2 else "Cac"
                elif atom.elem == "N":
                    free_os = len([1 for n in atom.neighbors if n.elem == "O" and len(n.neighbors) == 1])
                    atom.atom_type = "Ntr" if free_os >= 2 else "Npl"
            elif len(atom.neighbors) == 2:
                if atom.elem == "C":
                    atom.atom_type = "C1"
                elif atom.elem == "N":
                    atom.atom_type = "N2"
                elif atom.elem == "O":
                    atom.atom_type = "O3"
                elif atom.elem == "S":
                    atom.atom_type = "S3"

        # Valence 1 atoms only
        for atom in self.atoms:
            if len(atom.neighbors) != 1:
                continue
            n = atom.neighbors[0]
            if atom.elem == "C":
                atom.atom_type = "C1-" if n.elem == "O" and len(n.neighbors) == 1 else "C1"
            elif atom.elem == "N":
                atom.atom_type = "N1"
            elif atom.elem == "O":
                if n.atom_type in ["Ntr", "N1+", "Cac"]:
                    atom.atom_type = "O2-"
                elif n.atom_type in ["Pac", "N3+", "Sac", "Son", "Sxd"]:
                    atom.atom_type = "O3-"
                elif n.elem == "C":
                    atom.atom_type = "O1+" if len(n.neighbors) == 1 else "O2"
            elif atom.elem == "S":
                atom.atom_type = "S2"
            elif atom.elem == "H":
                atom.atom_type = "HC" if n.elem == "C" else "H"
            elif atom.elem == "D":
                atom.atom_type = "DC" if n.elem == "C" else "D"

        # Check charge of nitrogens
        for atom in self.atoms:
            if atom.elem == "N" and atom.atom_type != "N3+":
                if all([n.atom_type in ["C3", "H", "D"] and n.rings == 0 for n in atom.neighbors]):
                    atom.atom_type = "N3+"

            if atom.atom_type == "C2":
                npls = [a for a in atom.neighbors if a.atom_type in ["Npl", "Ng+"]]
                if len(npls) >= 2 and all([n.rings == 0 for n in npls]):
                    for n in npls:
                        n.atom_type = "Ng+"

    def build_grid(self):
        self.grid.clear()
        for atom in self.atoms:
            self.grid[(int(atom.pos[0]), int(atom.pos[1]), int(atom.pos[2]))].append(atom)

    def __init__(self, name, atom_records):
        """ Initialize a PDBModel object using a PDB-file or PDB-id. """

        self.name = name
        self.atoms = atom_records
        self.constraints = []
        self.grid = defaultdict(list)

        self.build_grid()
        self.build_cov_bonds()
        self.build_ring_counts()
        self.build_idatm()

    @staticmethod
    def hydrogen_bond_energy(aa, acceptor, hydrogen, donor):
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
        d0 = 8.0
        r0 = 2.8
        psi0 = 1.911  # 109.5degs
        r = acceptor.distance(donor)
        theta = angle(donor.pos, hydrogen.pos, acceptor.pos)
        psi = angle(hydrogen.pos, acceptor.pos, aa.pos)

        if r < 2.5 or r > 3.9:  # R > 3.2: #R > 3.9
            return 1000
        if theta < (3.141592/2):
            return 1000
        if donor.get_sp() == 3 and acceptor.get_sp() == 3 and psi-psi0 > (3.141592/2):
            return 1000
        if donor.get_sp() == 3 and acceptor.get_sp() == 2 and psi < (3.141592/2):
            return 1000

        ratio = r0/r
        ratio_pow2 = ratio*ratio
        ratio_pow4 = ratio_pow2*ratio_pow2
        ratio_pow8 = ratio_pow4*ratio_pow4
        ratio_pow10 = ratio_pow8*ratio_pow2
        ratio_pow12 = ratio_pow8*ratio_pow4

        energy_dist = d0 * (5*ratio_pow12 - 6*ratio_pow10)
        energy_angl = math.cos(theta)**2

        if donor.get_sp() == 3 and acceptor.get_sp() == 3:
            return energy_dist * energy_angl * math.cos(psi-psi0)**2
        if donor.get_sp() == 3 and acceptor.get_sp() == 2:
            return energy_dist * energy_angl * math.cos(psi)**2
        if donor.get_sp() == 2 and acceptor.get_sp() == 3:
            return energy_dist * energy_angl * energy_angl
        if donor.get_sp() == 2 and acceptor.get_sp() == 2:
            # compute phi
            n_d = normalize(cross(sub(donor.neighbors[0].pos, donor.pos), sub(donor.neighbors[1].pos, donor.pos)))
            if len(aa.neighbors) >= 2:  # Not necessarily the case
                n_a = normalize(cross(sub(aa.neighbors[0].pos, aa.pos), sub(aa.neighbors[1].pos, aa.pos)))
            else:
                n_a = normalize(cross(sub(aa.pos, acceptor.pos), sub(hydrogen.pos, acceptor.pos)))
            phi = math.acos(dot(n_d, n_a))
            if phi > 3.141592:
                phi = 3.141592-phi

            return energy_dist * energy_angl * math.cos(max(phi, psi))**2

        raise RuntimeError("Hybridization combination not implemented: " +
                           str(donor.get_sp())+" + "+str(acceptor.get_sp()))

    def compute_hydrogen_bonds(self, threshold=-1.0):
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
            if not acceptor.is_acceptor():
                continue
            for donor in self.get_approx_nearby(acceptor.pos, 3.9):
                if not donor.is_donor():
                    continue
                if acceptor == donor:
                    continue
                if len(acceptor.neighbors) == 0:
                    continue
                aa = acceptor.neighbors[0]
                for hydrogen in donor.neighbors:
                    if not hydrogen.elem == "H":
                        continue
                    if hydrogen.distance(acceptor) > 2.5:
                        continue
                    try:
                        energy = self.hydrogen_bond_energy(aa, acceptor, hydrogen, donor)
                        if energy < threshold:
                            bonds.append((aa, acceptor, hydrogen, donor, energy))
                    except:
                        pass
        sorted(bonds, key=lambda x: min(x[1].id, x[2].id))
        return bonds

    def compute_disulphidebonds(self):
        """Compute disulphide bonds by checking each cysteine for covalent neighbors."""
        ret = []
        sulphurs = [a for a in self.atoms if a.resn == "CYS" and a.elem == "S"]
        for s1 in sulphurs:
            s2 = [n for n in s1.neighbors if n.elem == "S"]
            if s2 and s2[0].id < s1.id:
                ret.append((s1, s2[0]))

    def burial(self, atom, radius=8.0):
        """Computes the amount of atoms within a `radius` from `atom`"""
        neighborhood = [a for a in self.get_nearby(atom.pos, radius)]
        return len(neighborhood)

    def __repr__(self):
        return "PDBModel('"+self.name+"')"

    def check_hydrogens(self):
        hydrogens = [a for a in self.atoms if a.elem == "H"]
        if not hydrogens:
            print("Warning: No hydrogens present in "+self.name+". Consider running `reduce` first.")

    def clean_alt_conformations(self, altloc):
        """
        Check all alternate conformations and prune them depending on the `altloc` parameter.

        Args:
            altloc (string): Only atoms with this alternate conformation will be retained, unless it has the value 'all'
        """
        if altloc == "all":
            return
        alt_confs = [a for a in self.atoms if not (a.alt in [" ", altloc])]
        self.atoms = [a for a in self.atoms if a.alt in [" ", altloc]]
        for a in self.atoms:
            a.alt = " "  # remove alternate conformation indicator
        if alt_confs:
            print(self.name+": Pruning "+str(len(alt_confs))+" atoms in alternative conformations")
        if alt_confs and verbose:
            for a in alt_confs:
                print("> removed alt conf "+str(a))

    def clean_waters(self):
        """Removes all HOH waters and prints a message"""

        waters = [a for a in self.atoms if a.resn == "HOH"]
        self.atoms = [a for a in self.atoms if a.resn != "HOH"]
        if waters:
            print(self.name+": Pruning "+str(len(waters))+" water atoms")
        if waters and verbose:
            for a in waters:
                print("> removed water "+str(a))

    def rebuild_atom_ids(self):
        """Renumber atom ids so they are consecutive and start at 1"""
        for i, a in enumerate(self.atoms):
            a.id = i+1

    def clean_residue_sequence(self):
        currentresid = self.atoms[0].resi
        currentchain = self.atoms[0].chain
        reslist = []  # temporarily stores atoms of one residue
        inconsistency_flag = 0
        for atom in self.atoms:
            if atom.hetatm:  # Only need continuous sequence for residues
                continue
            if atom.chain != currentchain:  # Reached a chain border, don't need to check for continous IDs
                currentchain = atom.chain  # Set to new chain
                currentresid = atom.resi  # Set to new first resi
                del reslist[:]  # clear list
                reslist.append((atom.name, atom.alt))
                if atom.alt != " ":  
                    reslist.append((atom.name, " "))  # Fill list with the "default" to correctly identify new residues
                if atom.alt == " ":
                    reslist.append((atom.name, "A"))  # Fill list with the "default" to correctly identify new residues
                continue
            if atom.resi != currentresid or ((atom.name, atom.alt) in reslist):  # New residue or wrong numbering
                # print atom.resi, atom.alt
                if atom.resi == currentresid + 1:  # desired ID in contiguous sequence, go on
                    currentresid = atom.resi  # move current residue one up
                    del reslist[:]  # clear list
                    reslist.append((atom.name, atom.alt))
                    if atom.alt != " ":  
                        reslist.append((atom.name, " "))
                    if atom.alt == " ":
                        reslist.append((atom.name, "A"))
                    continue
                if atom.resi > currentresid + 1:  # This is a gap
                    print("Chain "+str(currentchain)+": residue gap between "+str(currentresid)+" and "+str(atom.resi))
                    currentresid = atom.resi  # move current residue across gap (user has to ensure correct sequence)
                    del reslist[:]  # clear list
                    reslist.append((atom.name, atom.alt))
                    if atom.alt != " ":  
                        reslist.append((atom.name, " "))
                    if atom.alt == " ":
                        reslist.append((atom.name, "A"))
                    continue
                else:  # weird numbering (e.g. 52 and 52A)
                    if (atom.name, atom.alt) in reslist:  # New residue (even if same wrong number, or same name)
                        if inconsistency_flag == 0:
                            print("Inconsistent residue name at " +
                                  str(atom.resi)+atom.resn+" in chain "+str(currentchain))
                            inconsistency_flag = 1
                        currentresid += 1
                        atom.resi = currentresid
                        del reslist[:]  # clear list
                        reslist.append((atom.name, atom.alt))
                        if atom.alt != " ":  
                            reslist.append((atom.name, " "))
                        if atom.alt == " ":
                            reslist.append((atom.name, "A"))
                    else:  # wrong numbering, in same residue
                        atom.resi = currentresid  # adopt current ID
                        reslist.append((atom.name, atom.alt))
                        if atom.alt != " ":  
                            reslist.append((atom.name, " "))
                        if atom.alt == " ":
                            reslist.append((atom.name, "A"))
            else:      
                # Everything regular, just add name to list
                reslist.append(atom.name)  # keep track of atom names (each atom only once)
        
    def check_collisions(self):
        serious_collisions = []
        collisions = []
        for atom1 in self.atoms:
            for atom2 in self.get_approx_nearby(atom1.pos, atom1.vdw_radius()+2.0):  # 2.0 is the largest vdw radius for atom2
                if atom1.id <= atom2.id:
                    continue
                if abs(atom1.resi-atom2.resi) <= 1 and atom1.chain == atom2.chain:
                    continue
                dist = atom1.distance(atom2)
                if dist < 0.7*(atom1.covradius()+atom2.covradius()):
                    serious_collisions.append((atom1, atom2, dist))
                    collisions.append((atom1, atom2, dist))
                elif dist < 0.6*(atom1.vdw_radius()+atom2.vdw_radius()):
                    collisions.append((atom1, atom2, dist))
        if collisions or serious_collisions:
            print("%s: %d collisions, %d very serious" % (self.name, len(collisions), len(serious_collisions)))
            if verbose:
                for a1, a2, d in collisions:
                    print("> distance from %s to %s : %.2f" % (str(a1), str(a2), d))

        return collisions

    def check_covalent_bonds(self):
        """ Check if number of covalent neighbors make sense and print warnings if they don't. """
        odd_atoms = []
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "C" and not (len(a.neighbors) in [1, 2, 3, 4])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "N" and not (len(a.neighbors) in [1, 2, 3, 4])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "O" and not (len(a.neighbors) in [1, 2])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "H" and not (len(a.neighbors) in [1])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "S" and not (len(a.neighbors) in [2, 3, 4, 5, 6])]
        odd_atoms = odd_atoms+[a for a in self.atoms if a.elem == "P" and not (len(a.neighbors) in [2, 3, 4, 5])]
        if odd_atoms:
            print("%s: %d atoms with irregular number of covalent neighbors" % (self.name, len(odd_atoms)))
        if odd_atoms and verbose:
            for a in odd_atoms:
                print("> %s has %d covalent neighbors" % (str(a), len(a.neighbors)))

        return odd_atoms

    def pdb_str(self):
        """ Return a string representing this PDB model """
        ret = ""
        for atom in self.atoms:
            ret += atom.pdb_str()

        # Write CONECT records based on bond profile
        for atom in self.atoms:
            ret += atom.get_conect()

        return ret
      

    def write_pdb(self, fname):
        """ Write the model to a PDB file """
        f = open(fname, 'w')
        f.write(self.pdb_str())
        f.close()
    
    def get_atom(self, chain, resi, name):
        ret = [a for a in self.atoms if a.chain == chain and a.resi == resi and a.name == name]
        if ret:
            return ret[0]
        return None

    def getAtomsInResi(self, resi):
        ret = []
        for atom in self.atoms:
            if atom.resi==resi:
                ret.append(atom)
        return ret

    def getResidues(self):
        '''
        Return a sorted list of all residue numbers in this structure
        '''
        ret = set()
        for atom in self.atoms:
            ret.add(atom.resi)
        return sorted(list(ret))

    def getResidueIDsandNames(self):
        '''
        Return a sorted list of all residue numbers and names in this structure
        '''
        ret = set()
        for atom in self.atoms:
            ret.add((atom.resi,atom.resn))
        return sorted(list(ret))

class PDBFile:

    def pdbDownload(self,pdb_id):
        hostname="ftp.wwpdb.org"
        directory="/pub/pdb/data/structures/all/pdb/"
        prefix="pdb"
        suffix=".ent.gz"

        import os, sys, ftplib, shutil, gzip

        # Log into server
        #print "Downloading %s from %s ..." % (pdb_id, hostname)
        ftp = ftplib.FTP()
        ftp.connect(hostname)
        ftp.login()

        # Download all files in file_list
        to_get = "%s/%s%s%s" % (directory,prefix,pdb_id.lower(),suffix)
        to_write = "%s%s" % (pdb_id,suffix)
        final_name = "%s.pdb" % to_write[:to_write.index(".")]
        try:
            ftp.retrbinary("RETR %s" % to_get,open(to_write,"wb").write)
            f = gzip.open(to_write,'r')
            g = open(final_name,'w')
            g.writelines(f.readlines())
            f.close()
            g.close()
            os.remove(to_write)
        except ftplib.error_perm:
            os.remove(to_write)
            print "ERROR! %s could not be retrieved from PDB!" % to_get
            ftp.quit()
            return None

        # Log out
        ftp.quit()
        return final_name
        
    def __init__(self, pdb_file):
        """ Initialize a multi-model object using a PDB file name. """

        if len(pdb_file)==4: 
            pdb_file = self.pdbDownload(pdb)
            
        if pdb_file.endswith(".gz"): 
            import gzip
            f = gzip.open(pdb_file, 'rb')
        else:
            f = open(pdb_file, 'r')

        models = []
        atoms = []
        for line in f.readlines():
            if line[0:5] == "MODEL":
                if atoms:
                    raise RuntimeError("A MODEL was not stopped (ENDMDL) before a new one started")
            elif line[0:6] == "ENDMDL":
                if not atoms:
                    raise RuntimeError("Reached ENDMDL without reading any atoms")
                models.append(atoms)
                atoms = []
            elif line[0:4] == "ATOM" or line[0:6] == "HETATM":
                atom = Atom(line)
                atoms.append(atom)

        if atoms:
            models.append(atoms)

        f.close()
        base_name = pdb_file[0:pdb_file.rfind(".")]  # finds last ".", either .pdb or .ent
        if len(models) == 1:
            self.models = [PDBModel(base_name, models[0])]
        else:
            self.models = [PDBModel(base_name+"_"+str(i), m) for i, m in enumerate(models)]

    def rebuild_covalent_bonds(self):
        for m in self.models:
            m.rebuild_cov_bonds()

    def rebuild_atomids(self):
        for m in self.models:
            m.rebuild_atom_ids()

    def clean_dehydrohoh(self):
        for m in self.models:
            m.clean_dehydro_hoh()

    def clean_alt_conformations(self, altloc="all"):
        for m in self.models:
            m.clean_alt_conformations(altloc)

    def clean_waters(self):
        for m in self.models:
            m.clean_waters()
            m.build_grid()

    def clean_residue_sequences(self):
        for m in self.models:
            m.clean_residue_sequence()

    def check_hydrogens(self):
        for m in self.models:
            m.check_hydrogens()

    def check_collisions(self):
        for m in self.models:
            m.check_collisions()

    def check_covalent_bonds(self):
        for m in self.models:
            m.check_covalent_bonds()
            
    def getAtomsInResi(self, resi, model_number = 0):
        '''
        Return all atoms in residue
        '''
        m = self.models[model_number]
        return m.getAtomsInResi(resi)

    def getResidues(self, model_number = 0):
        '''
        Return a sorted list of all residue numbers in this structure
        '''
        m = self.models[model_number]
        return m.getResidues()
            
    def getResidueIDsandNames(self, model_number = 0):
        '''
        Return a sorted list of all residue numbers and names in this structure
        '''
        m = self.models[model_number]
        return m.getResidueIDsandNames()
