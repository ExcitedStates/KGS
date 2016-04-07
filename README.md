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
```
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release <path to>/source
$ make
```
If modifying files in the source simply type `make` in the build-dir to
recompile.  To use Intels compiler replace the third step with
```
$ CXX=icpc cmake -DCMAKE_BUILD_TYPE=Release <path to>/source
```

To compile with debugging symbols and optimizations disabled, open a terminal
and type
```
$ mkdir debug
$ cd debug
$ cmake -DCMAKE_BUILD_TYPE=Debug <path to>/source
$ make
```
This allows debuggers and analysers like gdb, lldb, or valgrind to give
meaningful source-code info.


## Usage 

TODO: Write when we've settled on cmd-line arguments
