
#ifndef QRTRANSPOSE_H
#define QRTRANSPOSE_H

#include <gsl/gsl_matrix.h>
#include "QR.h"


/**
 * Compute and return the QR decomposition of a given matrix' transpose. If
 * `M` is passed to the constructor, the `getQ` function will return the Q
 * matrix in the decomposition M^T = QR et similis for `getR`.
 *
 * Convenient for computation of QR-backed nullspace (see NullspaceQR)
 */
class TransposeQR{
 public:
  TransposeQR(gsl_matrix* M);

  void updateFromMatrix();

  gsl_matrix* getMatrix() const;

  /** Return the Q-matrix from the decomposition M^T = QR */
  gsl_matrix* getQ() const;

  /** Return the R-matrix from the decomposition M^T = QR */
  gsl_matrix* getR() const;

 private:
  const QR* m_qr;
  const gsl_matrix* m_origMatrix;
};

#endif
