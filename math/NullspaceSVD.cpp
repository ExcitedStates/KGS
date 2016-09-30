
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>

#include "NullspaceSVD.h"
#include "Logger.h"
#include "gsl_helpers.h"
#include "GSLQR.h"
#include "MKLQR.h"

double SINGVAL_TOL = 1.0e-12; //0.000000000001; // only generic 10^-12
double RDIAVAL_TOL = 1.0e-6;

using namespace std;

NullspaceSVD::NullspaceSVD(SVD * svd) :
  Nullspace(svd->matrix),
  m_svd(svd)
{
}

NullspaceSVD::~NullspaceSVD()
{
}

void NullspaceSVD::UpdateFromMatrix()
{
  m_svd->UpdateFromMatrix();

  //Compute nullspacesize
  double maxSingularValue = gsl_vector_get(m_svd->S, 0);

  //Case with m < n and all singular values non-zero
  m_nullspaceSize = std::max((int) (m_svd->V->size2 - m_svd->matrix->size1), 0);

  for (int i = 0; i < m_svd->S->size; ++i) {
    if (gsl_vector_get(m_svd->S, i) / maxSingularValue < SINGVAL_TOL) {
      m_nullspaceSize = m_svd->V->size2 - i;
      break;
    }
  }

  //TODO: If an existing basis of proper size is allocated we might not need to reallocate here
  if (m_nullspaceBasis)
    gsl_matrix_free(m_nullspaceBasis);

  if (m_nullspaceSize > 0) {
    gsl_matrix_view nullspaceBasis_view = gsl_matrix_submatrix(m_svd->V,
                                                               0,                               //Row
                                                               m_svd->V->size2 - m_nullspaceSize, //Col
                                                               m_svd->V->size2,                 //Height
                                                               m_nullspaceSize);                  //Width

    m_nullspaceBasis = gsl_matrix_calloc(m_svd->V->size2, m_nullspaceSize);
    gsl_matrix_memcpy(m_nullspaceBasis, &nullspaceBasis_view.matrix);
  }
  else {
    m_nullspaceBasis = gsl_matrix_calloc(m_svd->V->size2, 1);//1-dim vector with zeros as entries
  }
}


void Nullspace::WriteMatricesToFiles(
    const std::string& jac_file,
    const std::string& null_file,
    const std::string& sval_file) const
{
  gsl_matrix_outtofile(m_svd->matrix, jac_file);
  gsl_matrix_outtofile(m_svd->V, null_file);
  gsl_vector_outtofile(m_svd->S, sval_file);
}

gsl_matrix *Nullspace::getBasis() const {
  return m_nullspaceBasis;
}

SVD *Nullspace::getSVD() const {
  return m_svd;
}

Nullspace* Nullspace::createNullspace(gsl_matrix* M)
{
  return new SVDNullspace(M);
}
