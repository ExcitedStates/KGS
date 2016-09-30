
#ifndef KGS_NULLSPACEQR_H
#define KGS_NULLSPACEQR_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "math/QR.h"
#include "math/Nullspace.h"

/**
 * An implementation of Nullspace backed by a QR decomposition
 */
class NullspaceQR: public Nullspace {
 public:
  /** Will construct a nullspace using a QR (transpose) decomposition */
  NullspaceQR(QRTranspose* qr);

  ~NullspaceQR();

  /** Update the Nullspace (and underlying SVD) to reflect an updated state of the matrix */
  void UpdateFromMatrix();

private:
  QRTranspose* m_qr;                  ///< SVD underlying this nullspace

  friend class Configuration;
};


#endif //KGS_nullptrSPACE_H
