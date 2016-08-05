
#ifndef KGS_nullptrSPACE_H
#define KGS_nullptrSPACE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <string>

#include "math/SVD.h"

/**
 * Computes, stores, and maintains the nullspace of a gsl_matrix.
 * Supports projections of any gradient represented as a gsl_vector onto the null-space
 * and also analyzing and updating rigidity. Note that the Nullspace object only
 * accurately reflects the matrix (SVD) if the UpdateFromMatrix function has been called.
 * This is not done on construction.
 *
 * An SVD-pointer is passed to the Nullspace on construction. This pointer is kept through
 * the life-time of the Nullspace object (but ownership is not assumed) and if the contents
 * of the matrix is changed, the Nullspace can be update to reflect the change using the
 * abstract UpdateFromMatrix function.
 */
class Nullspace {
 public:
  Nullspace(SVD * svd);
  ~Nullspace();

  /** Print the nullspace to standard output */
  void print() const;

  /** Projects a vector on the nullspace */
  void ProjectOnNullSpace (gsl_vector *to_project, gsl_vector *after_project) const;

  /** Analyzes which dihedrals and hydrogen bonds are rigidified by constraints */
  void RigidityAnalysis(gsl_matrix* HBondJacobian);

  /** Update the Nullspace (and underlying SVD) to reflect an updated state of the matrix */
  void UpdateFromMatrix();

  /** Return the nullspace size */
  int NullspaceSize() const { return m_nullspaceSize; }

  /** Return number of cycle DOFs */
  int getNumDOFs() const { return n; }

  /** Return number of dihedral DOFs that were rigidified in the last call to RigidityAnalysis. */
  int NumRigidDihedrals() const { return numRigidDihedrals; }

  /** Return number of h-bond dihedrals that were rigidified in the last call to RigidityAnalysis. */
  int NumRigidHBonds() const { return numRigidHBonds; }

  /** Return the underlying matrix. */
  gsl_matrix* Matrix() const{ return svd->matrix; }

  /** Return the basis of the nullspace as columns of a matrix */
  gsl_matrix *getBasis() const;

  /** Return the SVD of the nullspace as columns of a matrix */
  SVD *getSVD() const;

  /**
   * Returns true iff the angle specified by the argument is rigidified.
   * This result is only accurate if UpdateFromMatrix and RigidityAnalysis have both
   * been called.
   */
  bool IsAngleRigid(int angle_id){ return fabs(gsl_vector_get(rigidAngles, angle_id)-1.0)<0.001; }

  void WriteMatricesToFiles(const std::string& jac_file,
                            const std::string& null_file,
                            const std::string& sval_file) const;



private:
  SVD* svd;                    ///< SVD underlying this nullspace
  int m_nullspaceSize;         ///< Size of nullspace (rank of jacobian)
  int m, n;                    ///< Dimensions of underlying matrix (jacobian)

  gsl_matrix* m_nullspaceBasis;///< Basis of the nullspace

  gsl_vector* rigidAngles;     ///< Binary vector indicating which m_dofs are rigid
  gsl_vector* rigidHBonds;     ///< Binary vector indicating which h-bonds are rigid

  int numCoordinatedDihedrals; ///< Coordinated dihedrals
  int numRigidDihedrals;       ///< Rigid dihedrals
  int numRigidHBonds;          ///< Rigid hydrogen bonds



  /// These values have to be chosen according to the numerical analysis
//  static constexpr double SINGVAL_TOL = 1.0e-12; //0.000000000001; // only generic 10^-12
  static constexpr double RIGID_TOL =   1.0e-10; //0.0000000001; //depends on molecule, but 10^-10 seems a good fit!

  friend class Configuration;
};


#endif //KGS_nullptrSPACE_H
