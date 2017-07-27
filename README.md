# Kino-geometric sampling (KGS) source

KGS is an engine that generates conformational perturbations in biomolecules 
(protein, RNA, or ligands alike) by maintaining user-specified non-local 
constraints. Specifically, collisions are avoided while hydrogen-bonds and 
disulphide bonds are maintained. A set of applications are built on this 
engine:

* [kgs_explore](https://github.com/ExcitedStates/KGS/wiki/kgs_explore/): Generating conformational ensemble from a seed structure
* [kgs_transition](https://github.com/ExcitedStates/KGS/wiki/kgs_transition/): Directed sampling from initial structure to target structure
* [kgs_rigidity](https://github.com/ExcitedStates/KGS/wiki/kgs_rigidity/): Rigidity analysis of structure file
* kgs_deerprep: In progress
* [kgs_relative_transition](https://github.com/ExcitedStates/KGS/wiki/kgs_relative_transition/): Driving conformations.


## Dependencies

KGS depends as a minimum on the GNU Scientific Library (GSL) only. GSL is
available through most standard package managers. For Mac, download
[homebrew](http://brew.sh/) and run `brew install gsl`. 

If the [Intel MKL libraries](http://software.intel.com/mkl) can be found on the system then faster and
automatically parallelized SVD computations will be enabled as well. 


## Compiling

To compile KGS from source go to a terminal and type
```bash
git clone https://github.com/ExcitedStates/KGS.git
mkdir KGS/build
cd KGS/build
cmake -DCMAKE_BUILD_TYPE=Release ../src
make -j
```
After modifying files in `src` simply type `make` in the build-dir to
recompile.  To use Intels compiler (v17.0.3 or later) replace the 
third step with
```bash
CXX=icpc cmake -DCMAKE_BUILD_TYPE=Release <repo path>/source
```

To compile with debugging symbols and optimizations disabled, open a terminal
and type
```bash
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug <repo path>/source
make
```
This allows debuggers and analysers like gdb, lldb, or valgrind to give
meaningful source-code info.


### Compiling on Stanfords Sherlock cluster

Sherlock has an old (4.4) GCC compiler and no MKL libraries by default. To 
compile on Sherlock start an `sdev` session and from an empty directory type
```bash
module load gcc 
module load intel
CXX=`which g++` cmake <repo path>/src
make -j 16
```

### Compiling on SSRLs unix systems at SLAC (including sdc's)

This system has an old (4.4) GCC compiler by default so the newer version 5 
should be enabled before running `cmake`. This is not necessary when subsequently
calling `make`. 
```bash
scl enable devtoolset-4 bash
cmake <repo path>/src
make -j 16
```

## Usage 

Before running, its recommended to put the executables and 
`kgs_prepare.py` script somewhere pointed to by the `PATH` environment 
variable. A number of options that affect e.g. motion-planning and 
collision detection can be passed, but a minimum working example is as 
follows
```bash
wget https://files.rcsb.org/download/2ERL.pdb
kgs_prepare.py 2ERL.pdb
kgs_explore --initial 2ERL.kgs.pdb
```
This should create an `output` directory and put 10 generated conformations 
in it (the default number of samples to generate). The prepared `.kgs.pdb`
file has the same coordinates as the original, but atom-ids have been 
renumbered, and REMARK records indicating hydrogen-bonds have been added. 
