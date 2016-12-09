#ifdef __INTEL_MKL
#include <mkl.h>

void my_mkl_free_buffers(){
	mkl_free_buffers();
}

#endif
