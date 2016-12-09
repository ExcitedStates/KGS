#include "CULASVD.h"

#include <algorithm>

#ifdef __GPU_CULA
#include <cula.hpp>


void CULASVD::printCulaError(culaStatus status)
{
  char buf[256];
  if(!status) return;
  culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
  log("dimitar") << endl << buf << endl;

  exit(1);
}

int CULASVD::meetsMinimumCulaRequirements()
{
  int cudaMinimumVersion = culaGetCudaMinimumVersion();
  int cudaRuntimeVersion = culaGetCudaRuntimeVersion();
  int cudaDriverVersion = culaGetCudaDriverVersion();
  int cublasMinimumVersion = culaGetCublasMinimumVersion();
  int cublasRuntimeVersion = culaGetCublasRuntimeVersion();

  if(cudaRuntimeVersion < cudaMinimumVersion) {
      log("dimitar")<<"InfoCula)\t CUDA runtime version is insufficient; version "<<cudaMinimumVersion<<" or greater is required."<<endl;
      return 0;
  }

  if(cudaDriverVersion < cudaMinimumVersion) {
      log("dimitar")<<"InfoCula)\t CUDA driver version is insufficient; version "<<cudaMinimumVersion<<" or greater is required."<<endl;
      return 0;
  }

  if(cublasRuntimeVersion < cublasMinimumVersion) {
      log("dimitar")<<"InfoCula)\t CUDA runtime version is insufficient; version "<<cublasMinimumVersion<<" or greater is required."<<endl;
      return 0;
  }

  return 1;
}

//void nr_svd (gsl_matrix* jac, gsl_vector* s_out, gsl_matrix* v_out) {
void CULASVD::UpdateFromMatrix()
{
	//const int m = jac->size1;
	//const int n = jac->size2;

  const char jobu  = 'A';    // 'N'o columns of U are computed and returned in array U
  const char jobvt = 'A';    // 'A'll n rows of VT are returned in array VT

  const int lda = m;         // "ld"=leading/physical dimension of array in mem
  const int ldu = m;
  const int ldvt = n;

	double *s = new double [std::min(m,n)];   // singular values of A, S(i) >= S(i+1), dimension=min(m,n)
	double *u = new double [m*m];
	double *vt = new double [n*n];

	// Start timer
  double start_time, end_time, svd_time;
	CTKTimer timer;
	start_time = timer.getTimeNow();

 	// Jac[row-major] = Jac_t[column-major] format
 	gsl_matrix* jac_t = gsl_matrix_alloc(n,m);
 	gsl_matrix_transpose_memcpy(jac_t,jac);
  // Assume lda=jac_t->tda=m (which should be)!
  double *a = jac_t->data;

  // Initialize CULA
  culaStatus status;
 	if(meetsMinimumCulaRequirements()) {
   	status = culaInitialize();
    log("dimitar") << "InfoCula)\t Initializing CULA with status... " << culaGetStatusString(status) << endl;
 		if(status != culaNoError) {
      log("dimitar") << "InfoCula)\t CULA Initialization error. Exit. " << endl;
 			exit(1);
 		}
  } else {
    log("dimitar") << "InfoCula)\t CULA Minimum requirements not met. Exit. " << endl;
    exit(1);
 	}

  if(!a || !s || !u || !vt) {
    log("dimitar") << "InfoCula)\t Host side allocation error." << endl;
    status = culaInsufficientMemory;
	  exit(1);
  }

  // Do SVD
  status = culaGesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt);

	// End timer
  end_time = timer.getTimeNow();
  svd_time = end_time - start_time;
  log("dimitar") << "InfoSVD:\t svd_time   = " << svd_time   << " seconds."<< endl;

  // Check CULA errors
  if(status != culaNoError) printCulaError(status);

	// Populate s_out with the S singular values
  log("dimitar") << "InfoSVD)\t s_out Vector (m=" << s_out->size << ")" << endl;
	//for (int i=0;i<m;++i)
	for (int i=0;i<s_out->size;++i)
		gsl_vector_set ( s_out, i, s[i] );
	if (s_out->size<std::min(m,n)) {
    log("dimitar") << "InfoSVD)\t row_dim (m) > col_dim (n) of Jacobian." << endl;
    log("dimitar") << "InfoSVD)\t Print remaining singular values: " << endl;
		for (int i=s_out->size;i<std::min(m,n);++i)
      log("dimitar") << "\t s[" << i << "]=" << s[i] << " ";
    log("dimitar") << endl;
	}

	// Populate v_out with the V matrix values (vt column-major=v row-major)
  log("dimitar") << "InfoSVD)\t v_out Matrix (rows=" << v_out->size1 << " columns=" << v_out->size2 << ")" << endl;
  for ( int i=0; i<n; ++i )
    for ( int j=0; j<n; ++j )
  		gsl_matrix_set (v_out, i, j, vt[i*n+j] );

 	// Free memory
 	gsl_matrix_free(jac_t);
 	delete[] s;
 	delete[] u;
	delete[] vt;

  // Shutdown
  culaShutdown();
  log("dimitar") << "InfoCula)\t Shutting down CULA with status... " << culaGetStatusString(status) << endl;

}


#else

void CULASVD::UpdateFromMatrix()
{
  throw "CULASVD::UpdateFromMatrix error! CULA not supported. Install CUDA+CULA and compile with -D__GPU_CULA";
}

#endif

