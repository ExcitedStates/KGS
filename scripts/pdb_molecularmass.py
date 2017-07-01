#!/usr/bin/python
import pdb_structure
import sys
import os.path

atomic_weight = {
    'H':    1.008, 'HE':   4.002, 'LI':   6.940, 'BE':   9.012, 'B':   10.810, 'C':   12.011, 'N':   14.007, 'O':   15.999,
    'F':   18.998, 'NE':  20.180, 'NA':  22.989, 'MG':  24.305, 'AL':  26.981, 'SI':  28.085, 'P':   30.973, 'S':   32.060,
    'CL':  35.450, 'AR':  39.948, 'K':   39.098, 'CA':  40.078, 'SC':  44.955, 'TI':  47.867, 'V':   50.942, 'CR':  51.996,
    'MN':  54.938, 'FE':  55.845, 'CO':  58.933, 'NI':  58.693, 'CU':  63.546, 'ZN':  65.380, 'GA':  69.723, 'GE':  72.630,
    'AS':  74.921, 'SE':  78.971, 'BR':  79.904, 'KR':  83.798, 'RB':  85.468, 'SR':  87.620, 'Y':   88.905, 'ZR':  91.224,
    'NB':  92.906, 'MO':  95.950, 'RU': 101.070, 'RH': 102.905, 'PD': 106.420, 'AG': 107.868, 'CD': 112.414, 'IN': 114.818,
    'SN': 118.710, 'SB': 121.760, 'TE': 127.600, 'I':  126.904, 'XE': 131.293, 'CS': 132.905, 'BA': 137.327, 'LA': 138.905,
    'CE': 140.116, 'PR': 140.907, 'ND': 144.242, 'SM': 150.360, 'EU': 151.964, 'GD': 157.250, 'TB': 158.925, 'DY': 162.500,
    'HO': 164.930, 'ER': 167.259, 'TM': 168.934, 'YB': 173.045, 'LU': 174.967, 'HF': 178.490, 'TA': 180.947, 'W' : 183.840,
    'RE': 186.207, 'OS': 190.230, 'IR': 192.217, 'PT': 195.084, 'AU': 196.966, 'HG': 200.592, 'TL': 204.380, 'PB': 207.200,
    'BI': 208.980, 'TH': 232.038, 'PA': 231.035, 'U':  238.028
}


def domainWeight(atoms):
    """
    Molecular weight in kDa of a list of atoms

    :param atoms:
    :return:
    """
    mol_weight = 0
    for a in atoms:
        if a.elem in atomic_weight:
            mol_weight+=atomic_weight[a.elem]
        else:
            print("Warning: Unknown element ignored: '%s'"%(a.elem))

    avogadro = 6.022e23
    print("Molecular weight: %.3f kDa = %.3E g"%(mol_weight/1000, mol_weight/avogadro))
    return mol_weight/1000


def domainCOM(atoms, massWeighted=True):
    """
    Returns domain center of mass / geometry if massWeighted=False
    :param atoms:
    :param massWeighted:
    :return:
    """
    x, y, z = 0.0, 0.0, 0.0
    totmass=0.0
    for a in atoms:
        if massWeighted:
            if a.elem in atomic_weight:
                m = atomic_weight[a.elem.upper()]
                x += a.x * m
                y += a.y * m
                z += a.z * m
                totmass += m
        else:
            x += a.x
            y += a.y
            z += a.z

    # print totmass
    if massWeighted:
        return x / totmass, y / totmass, z / totmass
    else:
        return x / len(atoms), y / len(atoms), z / len(atoms)

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print("Usage: "+sys.argv[0]+" <pdb-file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    struc = pdb_structure.PDBFile(pdb_file)

    name = os.path.basename(pdb_file).replace(".pdb","")

    mol_weight = 0
    for a in struc.models[0].atoms:
        if a.elem in atomic_weight:
            mol_weight += atomic_weight[a.elem]
        else:
            print("Warning: Unknown element ignored: '%s'"%(a.elem))
    avogadro = 6.022e23

    print("Molecular weight: %.3f kDa = %.3E g"%(mol_weight/1000, mol_weight/avogadro))

  
        

