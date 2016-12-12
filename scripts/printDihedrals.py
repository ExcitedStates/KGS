#!/usr/bin/python
import Bio.PDB
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                   help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                   const=sum, default=max,
				                      help='sum the integers (default: find the max)')

args = parser.parse_args()
print args.accumulate(args.integers)

for model in Bio.PDB.PDBParser().get_structure("1HMP", "1HMP.pdb") :
	for chain in model :
		poly = Bio.PDB.Polypeptide.Polypeptide(chain)
		print "Model %s Chain %s" % (str(model.id), str(chain.id)),
		print poly.get_phi_psi_list()
