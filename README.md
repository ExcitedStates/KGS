# Kino-geometric sampling (KGS) source

TODO: Write a short intro to the method

## Dependencies

KGS depends as a minimum on the GNU Scientific Library (GSL) only. GSL is
available through most standard package managers. For Mac, download
[homebrew](http://brew.sh/) and run `brew install gsl`. 

If the Intel MKL libraries can be found on the system then faster and
automatically parallelized SVD computations will be enabled as well. 


## Compiling

To compile KGS from source go to a terminal and type
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release <path to>/source
make
```
If modifying files in the source simply type `make` in the build-dir to
recompile.  To use Intels compiler (as of 2016 it doesn't support regular C++11 regex and therefore doesn't compile) replace the third step with
```bash
CXX=icpc cmake -DCMAKE_BUILD_TYPE=Release <path to>/source
```

To compile with debugging symbols and optimizations disabled, open a terminal
and type
```bash
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug <path to>/source
make
```
This allows debuggers and analysers like gdb, lldb, or valgrind to give
meaningful source-code info.


## Compiling on Sherlock at Stanford

Sherlock has an old (4.4) GCC compiler and no MKL libraries by default. To compile on Sherlock start an `sdev` session and from an empty directory type
```bash
module load gcc 
module load intel
CXX=`which g++` cmake <path to>/source
make -j 16
```

## Usage 

TODO: Write when we've settled on cmd-line arguments
