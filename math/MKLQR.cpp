#include "MKLQR.h"

#include <algorithm>

#ifdef __INTEL_MKL
#include <mkl_lapack.h>

void MKLQR::updateFromMatrix()
{
}

#else

void MKLQR::updateFromMatrix()
{
    throw "MKLQR::updateFromMatrix error! MKL not supported. Install MKL and compile with -D__INTEL_MKL";
}

#endif

