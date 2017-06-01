#!/usr/bin/python
import pdb_structure
import sys
import os.path
import math
import numpy as np

atomic_weight = {
  'H' :  1.008, 'He':  4.002, 'Li':  6.940, 'Be':  9.012,'B' : 10.810, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
  'F' : 18.998, 'Ne': 20.180, 'Na': 22.989, 'Mg': 24.305,'Al': 26.981, 'Si': 28.085, 'P' : 30.973, 'S' : 32.060,
  'Cl': 35.450, 'Ar': 39.948, 'K' : 39.098, 'Ca': 40.078,'Sc': 44.955, 'Ti': 47.867, 'V' : 50.942, 'Cr': 51.996,
  'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,'Cu': 63.546, 'Zn': 65.380, 'Ga': 69.723, 'Ge': 72.630,
  'As': 74.921, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,'Rb': 85.468, 'Sr': 87.620, 'Y' : 88.905, 'Zr': 91.224,
  'Nb': 92.906, 'Mo': 95.950, 'Ru':101.070, 'Rh':102.905,'Pd':106.420, 'Ag':107.868, 'Cd':112.414, 'In':114.818,
  'Sn':118.710, 'Sb':121.760, 'Te':127.600, 'I' :126.904,'Xe':131.293, 'Cs':132.905, 'Ba':137.327, 'La':138.905,
  'Ce':140.116, 'Pr':140.907, 'Nd':144.242, 'Sm':150.360,'Eu':151.964, 'Gd':157.250, 'Tb':158.925, 'Dy':162.500,
  'Ho':164.930, 'Er':167.259, 'Tm':168.934, 'Yb':173.045,'Lu':174.967, 'Hf':178.490, 'Ta':180.947, 'W' :183.840,
  'Re':186.207, 'Os':190.230, 'Ir':192.217, 'Pt':195.084,'Au':196.966, 'Hg':200.592, 'Tl':204.380, 'Pb':207.200,
  'Bi':208.980, 'Th':232.038, 'Pa':231.035, 'U' :238.028 
}


if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "Usage: "+sys.argv[0]+" <pdb-file>"
        sys.exit(1)

    pdb_file = sys.argv[1]
    struc = pdb_structure.PDBFile(pdb_file) 

    name = os.path.basename(pdb_file).replace(".pdb","")

    avgpos = np.array([0,0,0])
    mol_weight = 0.0
    for a in struc.models[0]:
        if a.elem in atomic_weight:
            mol_weight += atomic_weight[a.elem]
            avgpos += a.pos * atomic_weight[a.elem]
        else:
            print("Warning: Unknown element ignored: '%s'"%(a.elem)) 
    avgpos /= mol_weight

    rg = 0.0
    n = 0
    for a in struc.models[0]:
        if a.elem in atomic_weight:
            diff = avgpos - a.pos
            rg += diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]
            n += 1
    rg = math.sqrt(rg/n)
  
    print("Radius of gyration: %.3f Angstrom"%(rg))

  
        

