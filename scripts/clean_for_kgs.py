#!/usr/bin/python -tt

import sys
import string

from Bio.PDB import PDBParser,PDBIO,Select


def read_structure(nfilein):
	"""
	Read the pdb file
	PDB code is set to filename[:-4]
	"""
	parser=PDBParser(PERMISSIVE=True)
	pdbcode=nfilein[:-4]
	structure=parser.get_structure(pdbcode,nfilein)
	return structure

def write_structure(structure, nfileout, chainname):
	"""
	Write the clean structure in a pdb file
	Remove AltLoc by selecting A
	Remove water molecules
	Renumber atoms starting at 1
	(Optional select chain)
	"""

	io=PDBIO()
	io.set_structure(structure)
	
	# foundHetAtm = 0
	## Selecting chain, removing altlocs and water, renumbering
	class AtomSelect(Select):
		def accept_atom(self, atom):
			i=1
			# if(atom.get_record_type()=='HETATM'):
			# 	foundHetAtm = 1
			if not chainname or atom.get_parent().get_parent().get_id()==chainname: 
				if not atom.get_parent().get_id()[0]=='W':
					if not atom.is_disordered() or atom.get_altloc()=='A':
						atom.set_altloc(' ')
						atom.set_serial_number(i)
						i+=1
						return 1
					else:
						return 0
	# if (foundHetAtm == 1):
	# 	print "Found HetAtm entries, might have to adapt residue profiles."
	io.save(nfileout, AtomSelect())	
	
def main():
	"""
	Cleaning script to prepare for kgs input
	Input: input and output file names (chain id optional w/ -c)
	"""
	args=sys.argv[1:]
	if not args:
		print "Usage: [-c chain ] pdbin pdbout";
		sys.exit(1)
	chainname=''
	if args[0]=='-c':
		chainname=args[1]
		del args[0:2]
	if len(args)<2:
		print "Error: must specify input and output file names"
		sys.exit(1)
	(filein,fileout)=(args[0],args[1])

	structure=read_structure(filein)
	write_structure(structure, fileout, chainname)	


if __name__ == "__main__":
	main()

