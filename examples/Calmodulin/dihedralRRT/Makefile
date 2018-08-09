# Option 1: h-bonds, dihedral rrt
KGSARGS = --initial ../init_allH.kgs.pdb --workingDirectory ./ --samples 10 --metric dihedral --planner dihedralrrt --maxRotation 0.05 --preventClashes true --gradient 1

# Option 2: all constraints, dihedral rrt
#KGSARGS = --initial ../init_wHydro.kgs.pdb --workingDirectory ./ --samples 10 --metric dihedral --planner dihedralrrt --maxRotation 0.05 --preventClashes true --gradient 1

KGS = ../../../build/kgs_explore

# Assumes that KGSARGS has been set
output.txt:
	mkdir -p output
	$(KGS) $(KGSARGS) | tee $@

clean:
	rm -rf output output.txt kgs_planner.log kgs_report.log kgs_debug.log *.pdb
	rm -rf constraintPlots distancePlots energyPlots comboAnalysis
