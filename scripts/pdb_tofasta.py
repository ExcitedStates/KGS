#!/usr/bin/env python
import pdb_structure
import sys
import os.path



if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "Usage: "+sys.argv[0]+" <pdb-file>"
        sys.exit(1)
    pdbFile = sys.argv[1]
    
    struc = pdb_structure.PDBFile(pdbFile) 

    name = os.path.basename(pdbFile).replace(".pdb","")

    from collections import defaultdict
    chain_map = defaultdict(list)
    for a in struc.models[0]:
        chain_map[a.chain].append(a)

    protresnmap={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','ASX':'B','CYS':'C','GLU':'E','GLN':'Q','GLX':'Z','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    for c in chain_map:
        print(">"+name+"_"+c)
        resi_map = { a.resi: protresnmap[a.resn] for a in chain_map[c] }
        seq = []
        for r in range(min(resi_map.iterkeys()), max(resi_map.iterkeys())+1):
            seq.append( "-" if not r in resi_map else resi_map[r] )
        print("".join(seq))
