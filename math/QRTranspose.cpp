
#include "QRTranspose.h"
#include <iostream>

using namespace std;

QRTranspose::QRTranspose(gsl_matrix* M):
    m_origMatrix(M),
    m_qr(gsl_matrix_alloc(M->size2, M->size1))
{}

QRTranspose::~QRTranspose(){
  gsl_matrix_free(m_qr->getMatrix());
  delete m_qr;
}

void QRTranspose::updateFromMatrix()
{
  gsl_matrix_transpose_memcpy(m_qr->getMatrix(), m_origMatrix);
  m_qr->updateFromMatrix();
}

gsl_matrix* QRTranspose::getMatrix() const{
  return m_origMatrix;
}

gsl_matrix* QRTranspose::getQ() const{
  return m_qr->getQ();
}

gsl_matrix* QRTranspose::getR() const{
  return m_qr->getR();
}

