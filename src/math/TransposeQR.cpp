
#include "TransposeQR.h"
#include <iostream>

using namespace std;

TransposeQR::TransposeQR(gsl_matrix* M):
    m_origMatrix(M),
    m_qr(QR::createQR(gsl_matrix_alloc(M->size2, M->size1)))
{}


void TransposeQR::updateFromMatrix()
{
  gsl_matrix_transpose_memcpy(m_qr->getMatrix(), m_origMatrix);
  m_qr->updateFromMatrix();
}

gsl_matrix* TransposeQR::getMatrix() const{
  return m_origMatrix;
}

gsl_matrix* TransposeQR::getQ() const{
  return m_qr->getQ();
}

gsl_matrix* TransposeQR::getR() const{
  return m_qr->getR();
}

