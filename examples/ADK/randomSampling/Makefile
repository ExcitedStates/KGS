KGS = ../../../build/kgs_explore

#Dihedral RRT, collapsing rigid edges, using h-bonds and hydrophobics
KGSARGS = --initial ../init_withHydrophobics.kgs.pdb --collisionFactor 0.75 --seed 418 --saveData 2 --sampleRandom true --gradient 1 --maxRotation 0.05 --stepSize 0.1 --decreaseSteps 2 --workingDirectory ./ --preventClashes true --metric dihedral --planner dihedralrrt --samples 100 --collapseRigidEdges 1

# Assumes that KGSARGS has been set
all: output.txt

output.txt:
	mkdir -p output
	$(KGS) $(KGSARGS) | tee $@

memcheck:
	mkdir -p output
	rm -f valgrindOutput.txt
	valgrind --tool=memcheck --dsymutil=yes --leak-check=no --show-reachable=yes --num-callers=20 --track-fds=yes --track-origins=yes --log-file=valgrindOutput.txt $(KGS) $(KGSARGS)
	@echo "!! View valgrindOutput.txt to see output of memory check.!!"
	@echo "!! If you cant see line numbers copy the hex-code and paste !!"
	@echo "!! it in gdb. E.g.: break *0x1000AC502!!"

leakcheck:
	mkdir -p output
	rm -f valgrindOutput.txt
	valgrind --tool=memcheck --leak-check=yes --dsymutil=yes --show-reachable=yes --num-callers=20 --track-fds=yes --track-origins=yes --log-file=valgrindOutput.txt $(KGS) $(KGSARGS)
	@echo "!! View valgrindOutput.txt to see output of memory check !!"
	@echo "!! If you cant see line numbers copy the hex-code and paste !!"
	@echo "!! it in gdb. E.g.: break *0x1000AC502 !!"

clean:
	rm -rf output output.txt kgs_planner.log kgs_report.log kgs_debug.log
