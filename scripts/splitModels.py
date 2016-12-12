#!/usr/bin/python
import sys



def extractModels(pdbFile):
	"""
	Splits the multi-model PDB-file into several models
	"""
	model = -10;
	baseName = pdbFile[0:pdbFile.rfind(".")]
	outFile = None
	f = open(pdbFile, "r")
	for line in f:
		if line.startswith("MODEL"):
			model = int(line.strip()[line.strip().rfind(" ")+1:]);
			outName = baseName+"_"+str(model)+".pdb"
			outFile = open(outName, "w")
			print "Writing model "+str(model)+" to "+outName
			continue
		if model!=-10 and line.startswith("ENDMDL"):
			model = -10
			outFile.close()
			outFile = None
			continue
		if model!=-10 and line.startswith("ATOM"):
			outFile.write(line);
	f.close()


def main():
	"""
	Checks arguments and splits the multi-model PDB-file into several models
	"""

	if len(sys.argv)!=2:
		print "Usage: "+sys.argv[0]+" <pdb-file>"
		print "All models will be extracted from the specified pdb file (f.pdb) and written to f_m.pdb where m is each models number. Only ATOM records are printed."
		sys.exit(1)

	extractModels(sys.argv[1])

if __name__ == "__main__":
	main()
