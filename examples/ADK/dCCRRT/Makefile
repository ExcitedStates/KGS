KGS = ../../../build/kgs_transition
KGSARGS = --initial ../init.pdb --target ../target_aligned.pdb --hbondMethod user --hbondFile ../hBonds.txt --collisionFactor 0.75 --seed 418 --saveData 2 --sampleRandom false --gradient 3 --maxRotation 0.05 --stepSize 0.1 --decreaseSteps 2 --workingDirectory ./ --root -1 --preventClashes true --metric rmsd --planner dccrrt --samples 1000

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
	rm -rf output output.txt
