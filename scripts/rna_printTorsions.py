#!/usr/bin/python

from pdb_structure import *
import sys

def torsion(a1, a2, a3, a4):
	""" 
	Taken from 
	http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
	"""
	from numpy import cross, dot
	from numpy.linalg import norm
	import math
	b1 = (a2-a1)
	b2 = (a3-a2)
	b3 = (a4-a3)
	n1 = cross(b1,b2)/norm( cross(b1,b2) )
	n2 = cross(b2,b3)/norm( cross(b2,b3) )
	m1 = cross( n1, b2/norm(b2) )
	x = dot(n1, n2)
	y = dot(m1, n2)
	angle = -math.atan2(y,x)
	return angle



def printTorsions(pdb_file_name):
	pdb_file = PDBFile(pdb_file_name)

	first_res = 1000000
	last_res = -1000000
	for a in pdb_file.getAtoms():
		if a.resi<first_res: first_res = a.resi
		if a.resi>last_res: last_res = a.resi

	printFormat = "#%"+str(len(pdb_file_name)-1)+"s %4s %7s %7s %7s %7s %7s %7s %7s"
	print printFormat % ("ID","resi","alpha","beta","gamma","delta","epsilon","zeta","chi")

	for res in range(first_res, last_res+1):
		a1 = pdb_file.getAtom(res-1, "O3'")
		a2 = pdb_file.getAtom(res  , "P")
		a3 = pdb_file.getAtom(res  , "O5'")
		a4 = pdb_file.getAtom(res  , "C5'")
		a5 = pdb_file.getAtom(res  , "C4'")
		a6 = pdb_file.getAtom(res  , "C3'")
		a7 = pdb_file.getAtom(res  , "O3'")
		a8 = pdb_file.getAtom(res+1, "P")
		a9 = pdb_file.getAtom(res+1, "O5'")
		if not (a1 and a2 and a3 and a4 and a5 and a6 and a7 and a8 and a9):
			continue

		alpha = torsion(a1,a2,a3,a4)
		beta  = torsion(a2,a3,a4,a5)
		gamma = torsion(a3,a4,a5,a6)
		delta = torsion(a4,a5,a6,a7)
		epsil = torsion(a5,a6,a7,a8)
		zeta  = torsion(a6,a7,a8,a9)

		a10 = pdb_file.getAtom(res, "O4'")
		a11 = pdb_file.getAtom(res, "C1'")
		a12 = pdb_file.getAtom(res, "N9")
		a13 = pdb_file.getAtom(res, "C4")
		if not a12:
			a12 = pdb_file.getAtom(res, "N1")
			a13 = pdb_file.getAtom(res, "C2")

		chi = torsion(a10,a11,a12,a13)


		print "%s %4d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f" % (pdb_file_name,res,alpha,beta,gamma,delta,epsil,zeta,chi)


def main():

	if len(sys.argv)<=1:
		print "Usage: "+sys.argv[0]+" <pdb-file list>"
		sys.exit(1)

	for f in sys.argv[1:]:
		printTorsions(f)

if __name__ == "__main__":
	main()
