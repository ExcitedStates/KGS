#Option 1: with hydrogen bonds and hydrophobics (prepared with kgs_prepare)
#KGSARGS = --initial ../init_withHydrophobics.kgs.pdb --saveData 4 --workingDirectory ./

#Option 2: only hydrogen bonds (prepared with kgs_prepare)
#KGSARGS = --initial ../init.kgs.pdb --saveData 4 --workingDirectory ./

#Option 3: separate hydrogen bond file (for use with previous input)
#KGSARGS = --initial ../init.pdb --hbondMethod user --hbondFile ../hBonds.txt --saveData 4 --workingDirectory ./

#Option 4: also compute mutual information between two domains
KGSARGS = --initial ../init_withHydrophobics.kgs.pdb --saveData 4 --workingDirectory ./ --source "resi 30-59" --sink "resi 122-159"

include ../../Makefile_rigidity
