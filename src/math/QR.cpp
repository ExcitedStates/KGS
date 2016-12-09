#include <iostream>

#include <gsl/gsl_blas.h>

#include "QR.h"
#include "math/gsl_helpers.h"
#include "QRMKL.h"
#include "QRGSL.h"

using namespace std;

QR::QR(gsl_matrix* M):
    m_matrix(M),
    m(M->size1),
    n(M->size2),
    m_Q(gsl_matrix_alloc(m,m)),
    m_R(gsl_matrix_alloc(m,n))
{
}


QR::~QR(){
  gsl_matrix_free(m_Q);
  gsl_matrix_free(m_R);
}


void QR::print() const{
  std::cout<<"QR:"<<std::endl;
  std::cout<<"Q:"<<std::endl;
  gsl_matrix_cout(m_Q);
  std::cout<<"R:"<<std::endl;
  gsl_matrix_cout(m_R);
}


gsl_matrix* QR::getMatrix() const { return m_matrix; }
gsl_matrix* QR::getQ() const { return m_Q; }
gsl_matrix* QR::getR() const { return m_R; }


QR* QR::createQR(gsl_matrix* M)
{
#ifdef __INTEL_MKL
  return new QRMKL(M);
#else
  return new GSLQR(M);
#endif
}


