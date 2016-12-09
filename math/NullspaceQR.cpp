
#include "NullspaceQR.h"
#include "Logger.h"

double RDIAVAL_TOL = 1.0e-6;

using namespace std;

NullspaceQR::NullspaceQR(TransposeQR * qr) :
    Nullspace(qr->getMatrix()),
    m_qr(qr)
{
}

NullspaceQR::~NullspaceQR()
{
}

void NullspaceQR::updateFromMatrix()
{
  m_qr->updateFromMatrix();

  //Compute rank and nullspace size
  int rank = 0;
  for(int i=0;i<std::min(m,n);i++) {
    double val = gsl_matrix_get(m_qr->getR(), i, i);
    if (fabs(val) > RDIAVAL_TOL) rank++;
  }
  m_nullspaceSize = n-rank;

  //Free basis if already allocated
  if (m_nullspaceBasis)
    gsl_matrix_free(m_nullspaceBasis);

  //Extract nullspace from last columns of Q-matrix
  if (m_nullspaceSize > 0) {
    gsl_matrix_view nullspaceBasis_view = gsl_matrix_submatrix(m_qr->getQ(),
                                                               0,                                     //Row
                                                               n - m_nullspaceSize, //Col
                                                               n,                   //Height
                                                               m_nullspaceSize);                      //Width

    m_nullspaceBasis = gsl_matrix_calloc(n, m_nullspaceSize);
    gsl_matrix_memcpy(m_nullspaceBasis, &nullspaceBasis_view.matrix);

  }else {
    m_nullspaceBasis = gsl_matrix_calloc(m_matrix->size2, 1);//1-dim vector with zeros as entries
  }
}

