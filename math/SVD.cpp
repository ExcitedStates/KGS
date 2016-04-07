#include <iostream>

#include <gsl/gsl_blas.h>

#include "SVD.h"
#include "math/gsl_helpers.h"

using namespace std;

SVD::SVD(gsl_matrix* M):
    matrix(M),
    m(M->size1),
    n(M->size2),
    U(gsl_matrix_alloc(m,m)),
    S(gsl_vector_alloc(std::min(m,n))),
    V(gsl_matrix_alloc(n,n))
{
}


SVD::~SVD(){
  gsl_matrix_free(U);
  gsl_vector_free(S);
  gsl_matrix_free(V);
}

gsl_matrix* SVD::PseudoInverse() const{
  gsl_matrix* S_dag = gsl_matrix_calloc(n,m);
  const int sz = std::min(m,n);
  for(int i=0;i<sz;i++){
    double val = gsl_vector_get(S,i);
    if(val>0.0001) gsl_matrix_set(S_dag, i,i, 1/val);
  }

  gsl_matrix* prod1 = gsl_matrix_alloc(V->size1, S_dag->size2);
  gsl_matrix* prod2 = gsl_matrix_alloc(V->size1, U->size1);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, V, S_dag, 0.0, prod1);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, prod1, U, 0.0, prod2);

  gsl_matrix_free(prod1);
  gsl_matrix_free(S_dag);

  return prod2;
}

gsl_matrix* SVD::PseudoInverse(double lambda) const{
  gsl_matrix* S_dag = gsl_matrix_calloc(n,m);
  const int sz = std::min(m,n);
  for(int i=0;i<sz;i++){
    double val = gsl_vector_get(S,i);
    if(val>0.0001) gsl_matrix_set(S_dag, i,i, val/(val*val+lambda*lambda));
  }
  gsl_matrix* prod1 = gsl_matrix_alloc(V->size1, S_dag->size2);
  gsl_matrix* prod2 = gsl_matrix_alloc(V->size1, U->size1);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, V, S_dag, 0.0, prod1);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, prod1, U, 0.0, prod2);

  gsl_matrix_free(prod1);
  gsl_matrix_free(S_dag);

  return prod2;
}

//gsl_matrix* SVD::nullspace_projection()
//{
//  // Determine matrix rank, r
//  double epsilon = 10e-6;
//  int r = 0;
//  for(r=0;r<S->size; r++)
//    if( gsl_vector_get(S, r)<epsilon ) break;
//
//  // The rows from r to n of V^T represent the null-space basis transposed, N^T.
//  gsl_matrix_view N_t = gsl_matrix_submatrix(V_t, r,0, n-r, n);
//
//  // Projection matrix is defined as P = N N^T
//  gsl_matrix* P = gsl_matrix_alloc(n,n);
//  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &N_t.matrix, &N_t.matrix, 0.0, P);
//  return P;
//}

void SVD::print() const{
  std::cout<<"SVD:"<<std::endl;
  std::cout<<"U:"<<std::endl;
  gsl_matrix_cout(U);
  std::cout<<"S:"<<std::endl;
  gsl_vector_cout(S);
  std::cout<<"V:"<<std::endl;
  gsl_matrix_cout(V);
}
