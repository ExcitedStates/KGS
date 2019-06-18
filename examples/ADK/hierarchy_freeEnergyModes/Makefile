KGSARGS = --initial ../init_withHydrophobics.kgs.pdb --sampleFree true --samples 10 --saveData 1 --workingDirectory ./ 
KGS = /Users/StDoBudd/Documents/Forschung/Code/ExcitedStates/KGS/build/kgs_hierarchy 

all: output.txt 

output.txt: 
	 mkdir -p output
	 $(KGS) $(KGSARGS) | tee $@ 
clean:
	 rm -rf output output.txt hBonds_in.txt hBonds_out.txt hierarchy_data.txt