
#ifndef KGS_NULLSPACESVD_H
#define KGS_NULLSPACESVD_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "math/SVD.h"
#include "math/Nullspace.h"

/**
 * An implementation of Nullspace backed by SVD
 */
class NullspaceSVD: public Nullspace {
 public:
  /** Will construct a nullspace of `matrix` using the SVD decomposition */
  NullspaceSVD(SVD* svd);

  /** Update the Nullspace (and underlying SVD) to reflect an updated state of the matrix */
  void updateFromMatrix() override;

  /** Return the SVD of the nullspace as columns of a matrix */
  SVD *getSVD() const;

  void writeMatricesToFiles(
      const std::string& jac_file,
      const std::string& null_file,
      const std::string& sval_file) const;

private:
  SVD* m_svd;                  ///< SVD underlying this nullspace

  /// These values have to be chosen according to the numerical analysis
//  static constexpr double SINGVAL_TOL = 1.0e-12; //0.000000000001; // only generic 10^-12
  static constexpr double RIGID_TOL =   1.0e-10; //0.0000000001; //depends on molecule, but 10^-10 seems a good fit!

  friend class Configuration;
};


#endif //KGS_nullptrSPACE_H
