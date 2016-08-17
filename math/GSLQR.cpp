
#include <algorithm>
#include <iostream>
#include "GSLQR.h"

void GSLQR::updateFromMatrix()
{
  gsl_matrix* mCopy = gsl_matrix_alloc(m,n);
  gsl_matrix_memcpy(mCopy, m_matrix);

  gsl_vector* tau = gsl_vector_alloc(std::min(m,n));
  gsl_linalg_QR_decomp(mCopy, tau);
  gsl_linalg_QR_unpack(mCopy, tau, m_Q, m_R);

//  gsl_vector* tau = gsl_vector_alloc(std::min(m,n));
//  gsl_permutation* perm = gsl_permutation_alloc(n);
//  gsl_vector* norm = gsl_vector_alloc(n);
//  int* sign = new int(); *sign = 1;
//  gsl_linalg_QRPT_decomp2(m_matrix, m_Q, m_R, tau, perm, sign, norm );
//
//  std::cout<<"Permutation: ";
//  for(size_t i=0;i<n;i++) {
//    std::cout<<" "<<gsl_permutation_get(perm, i);
//  }
//  std::cout<<std::endl;
}
