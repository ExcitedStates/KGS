
#include <algorithm>
#include <iostream>
#include "GSLQR.h"

void GSLQR::updateFromMatrix()
{
  gsl_vector *tau = gsl_vector_alloc(std::min(m, n));
  gsl_permutation* perm = gsl_permutation_alloc(n);
  gsl_vector* norm = gsl_vector_alloc(n);
  int* sign = new int(); *sign = 1;
  gsl_linalg_QRPT_decomp2(m_matrix, m_Q, m_R, tau, perm, sign, norm );

  gsl_vector_free(tau);
  gsl_permutation_free(perm);
  gsl_vector_free(norm);
  free(sign);
}
