"""
Identifies CB atom pairs as potential locations for distance probes on given molecule
KGS must be installed for use
Creats a distribution_results folder with a result for each potential atom pair
Result given in the form of a histogram of the distribution of distances between the two atoms over time in the given ensemble
"""


from pdb_structure import *
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import sys


contact_radius = 5.0  # Max distance to neighboring residues
cutoff = 7  # Max number of neighboring residues to be considered exposed
total = 10  # Desired number of probe location options given as output


def getRigidAtoms():
    """Returns atoms in rigid bodies in PDB"""
    atoms = []
    fileName = pdbs[0]
    subprocess.call(["kgs_prepare.py", fileName])
    f = open(fileName[:-4] + ".kgs.pdb", "r")
    lineList = f.readlines()
    f.close()
    if connectivity and not "CONECT" in lineList[-1]:
        with open(connectivity) as fread:
            with open(fileName[:-4] + ".kgs.pdb", "a") as fwrite:
                for line in fread:
                    if "CONECT" in line:
                        fwrite.write(line)
    subprocess.call(["kgs_rigidity", "--initial", fileName[:-4] + ".kgs.pdb", "--saveData", "2", "--workingDirectory", "./"])
    with open(fileName[:-4] + ".kgs_RBs_1.txt") as f:
        for line in f:
            if not "NaN" in line and line != "\n":
                atoms.append(line[:-1])
    return atoms


def getNearbyRes(mol, atom):
    """Returns array of residues near given atom."""
    nearA = mol.get_nearby(atom, contact_radius)
    nearRes = []
    for a in nearA:
        if a.resi not in nearRes and a.resi != atom.resi:
            nearRes.append(a.resi)
    return nearRes


def getCBAtoms(mol):
    """Returns array of rigid CB atoms in given molecule."""
    CBs = []
    atoms = getRigidAtoms()
    for atom in mol.atoms:
        if atom.name == "CB" and str(atom.id) in atoms:
            CBs.append(atom)
    return CBs


def getSurfaceCBAtoms(mol):
    """Returns array of rigid surface CB atoms in given molecule."""
    result = []
    atoms = getCBAtoms(mol)
    for a in atoms:
        if len(getNearbyRes(mol, a)) < cutoff and a.chain == "A":
            result.append(a)
    return result


def getMolecules(pdbs):
    """Returns array of molecules from PDBs for distance data."""
    molecules = []
    i = 1
    for pdb in pdbs:
        print("opening " + str(i))
        tempmol = PDBFile(pdb).models[0]
        molecules.append(tempmol)
        i += 1
    return molecules


def getDistances(atom1, atom2):
    """Returns array of distances between two atoms."""
    distances = []
    for mol in molecules:
        a1 = getAtom(mol, atom1)
        a2 = getAtom(mol, atom2)
        d = a1.distance(a2)
        distances.append(d)
    return distances


def generateAtomKey(mol, atom):
    """Return key representation of atom."""
    key = str(mol) + ";" + str(atom.chain) + ";" + str(atom.resi) + ";" + str(atom.name)
    return key


def getResidueOfAtom(atom):
    """Return string form of residue of atom."""
    resstr = str(atom.chain) + " " + str(atom.resn) + " " + str(atom.resi)
    return resstr


def getAtom(mol, atom):
    """Gets atoms from cache (and adds to cache)."""
    key = generateAtomKey(mol, atom)
    if (key in atom_cache):
        result = atom_cache[key]
    else:
        result = mol.get_atom(atom.chain, atom.resi, atom.name)
        atom_cache[key] = result
    return result


if len(sys.argv) == 1:
    print("usage: " + sys.argv[0] + " <input>.pdb [options]")
    print("where <input> is every pdb in the ensemble")
    print("with option: --connectivity <filename>.pdb")
    print("    where <filename> is a pdb with CONECT lines to be used")
    sys.exit(-1)

connectivity = False
if "--connectivity" in sys.argv:
    connectivity = sys.argv[-1]
    sys.argv = sys.argv[:-2]
pdbs = sys.argv[1:len(sys.argv)]
molecules = getMolecules(pdbs)
atoms = getSurfaceCBAtoms(molecules[0])
sdevs = []  # Array for standard deviations of atom-pair distances
atompairs = []
atom_cache = {}
for an1 in range(len(atoms)):
    print(str(an1))
    for an2 in range(an1+1, len(atoms)):
        atompairs.append((an1, an2))
        distancesd = np.std(getDistances(atoms[an1], atoms[an2]))
        sdevs.append(distancesd)
sdevs_sorted = np.sort(sdevs)

directory = "distribution_results"
madeDir = False
i = 0
while not madeDir:
    if not os.path.exists(directory + "_" + str(i)):
        os.makedirs(directory + "_" + str(i))
        madeDir = True
    else:
        i += 1
num = total

min_distances = []
max_distances = []
best_atom_pairs = []
distances = []

for sdev in sdevs_sorted[-total:]:
    ind = sdevs.index(sdev)
    bestatoms = atompairs[ind]
    distanceDistribution = getDistances(atoms[bestatoms[0]], atoms[bestatoms[1]])
    best_atom_pairs.append(bestatoms)
    distances.append(distanceDistribution)
    min_distances.append(min(distanceDistribution))
    max_distances.append(max(distanceDistribution))

xmin = min(min_distances)
xmax = max(max_distances)

for atomnum in range(len(best_atom_pairs)):
    distanceDistribution = distances[atomnum]
    bestatoms = best_atom_pairs[atomnum]
    plt.hist(distanceDistribution, bins='auto', range=(xmin, xmax))
    plt.title("Distribution of distances between " + getResidueOfAtom(atoms[bestatoms[0]]) + " and " + getResidueOfAtom(atoms[bestatoms[1]]))
    plt.savefig(directory + "_" + str(i) + "/(" + str(num) + ") " + getResidueOfAtom(atoms[bestatoms[0]]) + " and " + getResidueOfAtom(atoms[bestatoms[1]]) + ".png", bbox_inches='tight')
    plt.close()
    num -= 1
