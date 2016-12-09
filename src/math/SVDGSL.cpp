
#include "SVDGSL.h"

void SVDGSL::UpdateFromMatrix()
{
  if (m<n) {
    //TODO: Memorize U_t, work tau and R for later calls to update
    gsl_matrix* U_t = gsl_matrix_alloc(n,m);
    gsl_matrix_transpose_memcpy(U_t,matrix);
    gsl_vector *work = gsl_vector_alloc(m);

    gsl_linalg_SV_decomp(U_t,U,S,work);

    // QR decomposition to extend U_t (the final V^T) to a full nxn basis
    gsl_vector *tau = gsl_vector_alloc(m);
    gsl_linalg_QR_decomp(U_t,tau);
    gsl_matrix *R = gsl_matrix_alloc(n,m);
    gsl_linalg_QR_unpack(U_t, tau, V, R);
    gsl_matrix_transpose(V);

    // free memory
    gsl_vector_free(work);
    gsl_matrix_free(U_t);
    gsl_matrix_free(R);
    gsl_vector_free(tau);
  }
  else {
    //TODO: Use QR unpack to get full U matrix.
    //The U matrix is not necessary for the null-space, but this code makes pseudo-inverses
    //fail when m>n.

    //gsl_matrix_free(U);
    //U = gsl_matrix_alloc(m,n);
    //gsl_matrix_memcpy(U,matrix);
    //gsl_vector * work = gsl_vector_alloc(n);
    //gsl_linalg_SV_decomp(U, V_t,S,work);
    //gsl_vector_free(work);
    gsl_vector * work = gsl_vector_alloc(n);
    gsl_matrix * U_tmp = gsl_matrix_alloc(m,n);
    gsl_matrix_memcpy(U_tmp,matrix);
    gsl_linalg_SV_decomp(U_tmp, V, S, work);
    gsl_vector_free(work);
    gsl_matrix_free(U_tmp);
  }
}
