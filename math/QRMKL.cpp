#include "QRMKL.h"
#include "gsl_helpers.h"

#include <algorithm>
#include <iostream>

#ifdef __INTEL_MKL
#include <mkl.h>

void QRMKL::updateFromMatrix()
{
  //Use m_R to store input to and then output from dgeqp3 temporarily
  gsl_matrix_memcpy(m_R, m_matrix);

  //Call LAPACKE_dgeqp3 on m_R (i.e. contents of m_matrix) and store the result in m_R
  int* jpvt = (int*)calloc(n, sizeof(int));
  int reflectors = std::min(m,n);
  double* tau = (double*)calloc(reflectors, sizeof(double));
  //Documentation for dgeqp3: https://software.intel.com/en-us/node/521004
  int status = LAPACKE_dgeqp3(LAPACK_ROW_MAJOR, m, n, m_R->data, n, jpvt, tau);
  if(status!=0) throw "QRMKL::updateFromMatrix error: Call to LAPACKE_dgeqp3 failed";

  //Use dormqr to multiply the implicitly encoded Q-matrix with the identity matrix to get the explicit Q
  gsl_matrix_set_identity(m_Q);
  //Documentation for dormqr: https://software.intel.com/en-us/node/521011
  status = LAPACKE_dormqr(LAPACK_ROW_MAJOR, 'L', 'N', m, m, reflectors, m_R->data, n, tau, m_Q->data, m);
  if(status!=0) throw "QRMKL::updateFromMatrix error: Call to LAPACKE_dormqr failed";

  //The upper triangle of m_R contains the actual R matrix and the lower triangle contains the
  //implicitly encoded Q matrix, so now simply reset the lower triangle to 0.0
  for(int i=1;i<m;i++) {
    std::fill_n(m_R->data + (i * n), std::min(i,n), 0.0);
  }

  free(tau);
  free(jpvt);
}

#else

void QRMKL::updateFromMatrix()
{
    throw "QRMKL::updateFromMatrix error! MKL not supported. Install MKL and compile with -D__INTEL_MKL";
}

#endif

