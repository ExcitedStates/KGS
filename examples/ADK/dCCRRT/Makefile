# joined list of h-bonds
#KGSARGS = --initial ../init.pdb --target ../target_aligned.pdb --hbondMethod user --hbondFile ../hBonds.txt --collisionFactor 0.75 --seed 418 --saveData 2 --sampleRandom false --gradient 3 --maxRotation 0.05 --stepSize 0.1 --decreaseSteps 2 --workingDirectory ./ --root -1 --alignIni true --preventClashes true --metric rmsd --planner dccrrt --samples 1000 --biasToTarget 0.4 --sampleReverse true

#separate list of h-bonds (in kgs.pdb files)
KGSARGS = --initial ../init.kgs.pdb --target ../target_aligned.kgs.pdb --collisionFactor 0.75 --seed 418 --saveData 2 --sampleRandom false --gradient 3 --maxRotation 0.05 --stepSize 0.1 --decreaseSteps 2 --workingDirectory ./ --root -1 --alignIni true --preventClashes true --metric rmsd --planner dccrrt --samples 100 --biasToTarget 0.4 --sampleReverse true

include ../../Makefile_transition
