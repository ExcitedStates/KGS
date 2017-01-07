#!/usr/bin/python

import sys
import os


def main():


	if len(sys.argv)<2:
		print "Usage: "+sys.argv[0]+"<input pdb file> to process [output filename, optional]"
		sys.exit(1)

	inputFile = sys.argv[1]
	outputFile1 = inputFile[:-4]+"_clean.pdb"
	# outputFile2 = inputFile[:-4]+"_Processed.pdb"
	outputFile2 = "init.pdb"
	hbFile = "hbonds_init.txt"
	if len(sys.argv)>2:
		outputFile2 = sys.argv[2]

	#Extract the first model
	print "Extracting the first model"
	getModel = "python /Users/StDoBudd/Documents/Forschung/Code/KGSrepo_SimTK/kgs/trunk/Scripts/extractModel.py "+inputFile+" 1 > "+outputFile2
	os.system(getModel)

	#Add hydrogen atoms
	print "Adding hydrogens using Reduce"
	runReduce = "/Users/StDoBudd/Documents/Forschung/Code/Reduce/reduce "+outputFile2+" > "+outputFile1
	os.system(runReduce)
	
	#Clean pdb for kgs (atom renumber, remove alt locs, remove waters)
	print "Cleaning file for KGS"
	runAtRenum = "python /Users/StDoBudd/Documents/Forschung/Code/KGSrepo_SimTK/kgs/trunk/Scripts/clean_for_kgs.py "+outputFile1+" "+outputFile2
	os.system(runAtRenum)
	
	#Run hbplus to identify hydrogen bonds
	# print "Identifying hydrogen bonds using hbplus"
	# runHBplus = "/Users/StDoBudd/Documents/Forschung/Code/hbplus/hbplus "+outputFile2
	# os.system(runHBplus)
	
	#Run hbfinder to identify hydrogen bonds
	print "Identifying hydrogen bonds using hbfinder"
	runHBfinder = "hbfinder -1 "+outputFile2+" > "+hbfile
	os.system(runHBfinder)

if __name__ == "__main__":
 	main()


