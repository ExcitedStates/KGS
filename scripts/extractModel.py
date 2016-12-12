#!/usr/bin/python
import sys



def extractModel(pdbFile, modelNum):
	"""
	Extracts the specified model from a multi-model PDB-file
	"""
	inModel = False;
	modelFile = False; #Does this file have separate models?
	f = open(pdbFile, "r")
	for line in f:
		if (line.startswith("MODEL") and modelFile == False):
			modelFile = True;
		if ( line.startswith("MODEL") and modelFile == True and line.strip().endswith(" "+modelNum) ):
			inModel = True;
			continue;
		if ( line.startswith("ATOM") and modelFile == False ): #Case without different models
			inModel = True;
		if inModel and ( line.startswith("ENDMDL") ):
			break;
		if inModel and ( line.startswith("ATOM") ): #or line.startswith("HETATM") or line.startswith("CONECT")):
			print line.strip();
	f.close()
	'''
	with open(pdbFile, "r") as f:
		for line in f:
			if line.startswith("MODEL") and line.strip().endswith(" "+modelNum):
				inModel = True;
				continue;
			if inModel and line.startswith("ENDMDL"):
				break;
			if inModel and line.startswith("ATOM"):
				print line.strip();
				'''


def main():
	"""
	Checks arguments and extracts the specified model from a multi-model PDB-file, default model: 1
	"""

	if len(sys.argv)<2:
		print "Usage: "+sys.argv[0]+" <pdb-file> <model-number = 1>"
		print "The specified model will be extracted from the pdb file and printed to std-out."
		print "ATOM, HETATM, and CONECT records are printed."
		sys.exit(1)

	modelNumber=str(1)
	if len(sys.argv)==3:
		modelNumer = sys.argv[2]
	extractModel(sys.argv[1],modelNumber)

if __name__ == "__main__":
	main()
