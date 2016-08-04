
#include <algorithm>
#include "GSLQR.h"

void GSLQR::updateFromMatrix()
{
  //gsl_matrix* matrixTrans = gsl_matrix_alloc(n,m);
  //gsl_matrix_transpose_memcpy(matrixTrans, m_matrix);
  gsl_matrix* mCopy = gsl_matrix_alloc(m,n);
  gsl_matrix_memcpy(mCopy, m_matrix);

  gsl_vector* tau = gsl_vector_alloc(std::min(m,n));
  gsl_linalg_QR_decomp(mCopy, tau);
  gsl_linalg_QR_unpack(mCopy, tau, m_Q, m_R);
  //Function: int gsl_linalg_QR_decomp (gsl_matrix * A, gsl_vector * tau)
  //Function: int gsl_linalg_QR_unpack (const gsl_matrix * QR, const gsl_vector * tau, gsl_matrix * Q, gsl_matrix * R)

}
