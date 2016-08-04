
#ifndef QR_H
#define QR_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

/**
 * Computes, stores, and maintains a QR decomposition of a gsl_matrix.
 *
 * A matrix-pointer is passed to the QR on construction. This pointer is kept through
 * the life-time of the QR object (but ownership is not assumed). If the contents
 * of the matrix is changed, the QR can be update to reflect the change using the
 * abstract UpdateFromMatrix function. This is not automatically done on construction.
 */
class QR{
 protected:
  const int m, n; ///< Dimensions of matrix

  /** Decomposes M into U*S*V^t.  */
  QR(gsl_matrix* M);
 public:


  /**
   * Constructs a QR object. If MKL is available it will be an MKLQR and if GSL
   * is available it will be GSLQR
   */
  static QR* createQR(gsl_matrix* M);

  virtual ~QR();

  /** Update Q adn R to reflect the decomposition of the matrix. */
  virtual void updateFromMatrix() = 0;

  /** Print the QR to standard out */
  void print() const;
 protected:
  gsl_matrix * const m_matrix;
  gsl_matrix * const m_Q;
  gsl_matrix * const m_R;

};

#endif
