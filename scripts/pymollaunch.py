#!/usr/bin/python -tt

import sys
import string
from os import sep

from glob import glob

import pymol

def launch(pdbcode):
	pymol.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
	pymol.finish_launching()

	pymol.cmd.set("cartoon_ring_mode", 3)
	pymol.cmd.set("cartoon_ring_finder", 1)
	pymol.cmd.set("cartoon_ladder_mode", 1)
	pymol.cmd.set("cartoon_nucleic_acid_mode", 4)
	pymol.cmd.set("cartoon_ring_transparency", 0.5)
	
	pymol.cmd.load(pdbcode+'_ok.pdb',pdbcode+'_init')
	pymol.cmd.color("hotpink",pdbcode+'_init')
	lst = glob("output/newpdb*.pdb")
	lst.sort()
#	for filen in lst: pymol.cmd.load(filen,"mov")
	for filen in lst:
		filenok=string.replace(filen,sep,'_')
		print filenok
		pymol.cmd.load(filen,filenok)
		pymol.cmd.color("slate",filenok)
		pymol.cmd.pair_fit(filenok,pdbcode+'_init')
	pymol.cmd.show("cartoon")
	pymol.cmd.hide("lines")

	pymol.cmd.center(pdbcode+'_init')
	pymol.cmd.save(pdbcode+'.pse')
	
def main():
	if len(sys.argv)==2:
		launch(sys.argv[1])		
	else:
		print "Usage: pdbcode"
	
if __name__ == '__main__':
	main()

