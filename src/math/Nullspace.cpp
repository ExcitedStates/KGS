
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "Nullspace.h"
#include "Logger.h"
#include "gsl_helpers.h"
//#include "QRGSL.h"
//#include "QRMKL.h"
#include "NullspaceSVD.h"
//#include "NullspaceQR.h"

using namespace std;

Nullspace::Nullspace(gsl_matrix* M) :
    m_matrix(M),
    m(M->size1),
    n(M->size2),
    numCoordinatedDihedrals( 0 ),
    numRigidDihedrals( 0 ),
    numRigidHBonds( 0 ),
    numRigidHydrophobicBonds( 0 ),
    m_nullspaceBasis(nullptr),
    m_rigidCovBonds(gsl_vector_alloc(n)),///for allocation, use maximum size of all n covalent edges
    m_rigidHBonds(gsl_vector_alloc(m)),///for allocation, use maximum size of all m constraints
    m_rigidDBonds(gsl_vector_alloc(m)),
    m_rigidHydrophobicBonds(gsl_vector_alloc(m))
{
}

Nullspace::~Nullspace ()
{
  if(m_nullspaceBasis)
    gsl_matrix_free(m_nullspaceBasis);
  gsl_vector_free(m_rigidCovBonds);
  gsl_vector_free(m_rigidHBonds);
  gsl_vector_free(m_rigidDBonds);
  gsl_vector_free(m_rigidHydrophobicBonds);
}


/*void Nullspace::performHydroRigidityAnalysis(gsl_matrix *HydrophobicBondJacobian)
{
  // First, check the dihedral angles for rigidity

 // gsl_vector_set_zero(m_rigidCovBonds);
  gsl_vector_set_zero(m_rigidHydrophobicBonds);

 numCoordinatedDihedrals = 0;
  numRigidDihedrals = 0;
  gsl_matrix* N = getBasis();

  for(int i=0; i<n; i++){
    bool moving=false;

    for(int j=0; j<N->size2; j++){
      double val = fabs( gsl_matrix_get(N,i,j) );
      if( val > RIGID_TOL ) {
        moving = true;
        break;
      }
    }

    // at least one entry was greater than threshold --> dihedral is coordinated
    if( moving ) {
      numCoordinatedDihedrals++;
    } else {
      gsl_vector_set(m_rigidCovBonds,i,1); /// binary list of rigid dihedrals, complement to coordinated version
      numRigidDihedrals++;
    }
  }

  log("constraints")<<"There are "<<numRigidDihedrals << " rigidified and " << numCoordinatedDihedrals << " coordinated dihedrals" << endl;

  // Now, check the hydrophobic m_bonds for rigidity
  int numHydrophobicBonds = HydrophobicBondJacobian->size1;
  gsl_matrix* HydrophobicBondNullspace = gsl_matrix_alloc(numHydrophobicBonds, std::max(m_nullspaceSize,1));
  gsl_vector* currentHydrophobicBondRow = gsl_vector_alloc(std::max(m_nullspaceSize, 1));

  ///Calculate the "In-Nullspace" Rotation of the hBonds
  //gsl_matrix_view m_nullspaceBasis = gsl_matrix_submatrix(svd->V,0,svd->V->size2-m_nullspaceSize,n,m_nullspaceSize); //Matrix, top-left element, num rows, num cols
  //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HBondJacobian, &m_nullspaceBasis.matrix, 0.0, hBondNullspace);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HydrophobicBondJacobian, m_nullspaceBasis, 0.0, HydrophobicBondNullspace);

  //Only for test analysis: write hBondNullspace=J_h*N to file
//	string outNullH="KINARI_Comparison/output/hBondnullSpace.txt";
//	gsl_matrix_outtofile(hBondNullspace, outNullH);

  double maxVal;
  double minVal;
  numRigidHydrophobicBonds = 0;

  for(int i=0; i<numHydrophobicBonds; i++){
    gsl_matrix_get_row(currentHydrophobicBondRow,HydrophobicBondNullspace,i);
    gsl_vector_minmax(currentHydrophobicBondRow, &minVal, &maxVal);

    if (minVal > -RIGID_TOL && maxVal < RIGID_TOL){
      gsl_vector_set(m_rigidHydrophobicBonds,i,1);
      numRigidHydrophobicBonds++;
    }
  }

  log("constraints")<<"There are "<<numRigidHydrophobicBonds<<" rigid out of "<<numHydrophobicBonds<<" hydrophobic bonds!"<<endl;

  gsl_vector_free(currentHydrophobicBondRow);
  gsl_matrix_free(HydrophobicBondNullspace);
}
*/
/** Analyzes which dihedrals and hydrogen bonds are rigidified by constraints */
void Nullspace::performRigidityAnalysis(gsl_matrix *HBondJacobian, gsl_matrix *DBondJacobian, gsl_matrix *HydrophobicBondJacobian)
{
  // First, check the dihedral angles for rigidity

  gsl_vector_set_zero(m_rigidCovBonds);
  gsl_vector_set_zero(m_rigidHBonds);
  gsl_vector_set_zero(m_rigidDBonds);
  gsl_vector_set_zero(m_rigidHydrophobicBonds);
  numCoordinatedDihedrals = 0;
  numRigidDihedrals = 0;
  gsl_matrix* N = getBasis();

  for(int i=0; i<n; i++){
    bool moving=false;

    for(int j=0; j<N->size2; j++){
      double val = fabs( gsl_matrix_get(N,i,j) );
      if( val > RIGID_TOL ) {
        moving = true;
        break;
      }
    }

    // at least one entry was greater than threshold --> dihedral is coordinated
    if( moving ) {
      numCoordinatedDihedrals++;
    } else {
      gsl_vector_set(m_rigidCovBonds,i,1); /// binary list of rigid dihedrals, complement to coordinated version
      numRigidDihedrals++;
    }
  }

  log("constraints")<<"There are "<<numRigidDihedrals << " rigidified and " << numCoordinatedDihedrals << " coordinated dihedrals" << endl;

  // Now, check all constraint bonds for their rigidity
  double maxVal;
  double minVal;
  /// HYDROGEN BOND CONSTRAINTS
  if(HBondJacobian) {/// there are hydrogen bonds
    int numHBonds = HBondJacobian->size1;
    gsl_matrix *hBondNullspace = gsl_matrix_alloc(numHBonds, std::max(m_nullspaceSize, 1));
    gsl_vector *currentHBondRow = gsl_vector_alloc(std::max(m_nullspaceSize, 1));

    ///Calculate the "In-Nullspace" Rotation of the hBonds
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, HBondJacobian, m_nullspaceBasis, 0.0, hBondNullspace);

    numRigidHBonds = 0;
    numRigidHydrophobicBonds = 0;
    for (int i = 0; i < numHBonds; i++) {
      gsl_matrix_get_row(currentHBondRow, hBondNullspace, i);
      gsl_vector_minmax(currentHBondRow, &minVal, &maxVal);

      if (minVal > -RIGID_TOL && maxVal < RIGID_TOL) {
        gsl_vector_set(m_rigidHBonds, i, 1);
        numRigidHBonds++;
      }
    }

    log("constraints") << "There are " << numRigidHBonds << " rigid out of " << numHBonds << " hydrogen bonds!" << endl;

    gsl_vector_free(currentHBondRow);
    gsl_matrix_free(hBondNullspace);
  }
    /// DISTANCE CONSTRAINTS
  if (DBondJacobian) {//there are distance bonds
    int numDBonds = (DBondJacobian->size1) / 3;
    gsl_matrix *dBondNullspace = gsl_matrix_alloc(DBondJacobian->size1, std::max(m_nullspaceSize, 1));
    gsl_vector *currentDBondRow = gsl_vector_alloc(std::max(m_nullspaceSize, 1));
    ///Calculate the "In-Nullspace" Rotation of the hydrophobic bonds
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, DBondJacobian, m_nullspaceBasis, 0.0, dBondNullspace);

    for (int i = 0; i < numDBonds; i++) {
      double minValThreeRows = 99999;
      double maxValThreeRows = 0;
      for (int j = 0; j < 3; j++) {//check all 3 remaining motions of the distance bond
        gsl_matrix_get_row(currentDBondRow, dBondNullspace, i * 3 + j);
        gsl_vector_minmax(currentDBondRow, &minVal, &maxVal);
        if (minVal < minValThreeRows)
          minValThreeRows = minVal;
        if (maxVal > maxValThreeRows)
          maxValThreeRows = maxVal;
      }

      if (minValThreeRows > -RIGID_TOL && maxValThreeRows < RIGID_TOL) {
        gsl_vector_set(m_rigidDBonds, i, 1);
        numRigidDBonds++;
      }
    }
    log("constraints") << "There are " << numRigidDBonds << " rigid out of " << numDBonds << " distance bonds!" << endl;
    gsl_vector_free(currentDBondRow);
    gsl_matrix_free(dBondNullspace);
  }
    /// HYDROPHOBIC CONSTRAINTS
  if (HydrophobicBondJacobian) {//there are hydrophobic bonds
    int numHydrophobicBonds = (HydrophobicBondJacobian->size1) / 5;
    gsl_matrix *hydrophobicBondNullspace = gsl_matrix_alloc(HydrophobicBondJacobian->size1,std::max(m_nullspaceSize, 1));
    gsl_vector *currentHydrophobicBondRow = gsl_vector_alloc(std::max(m_nullspaceSize, 1));
    ///Calculate the "In-Nullspace" Rotation of the hydrophobic bonds
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, HydrophobicBondJacobian, m_nullspaceBasis, 0.0, hydrophobicBondNullspace);

    for (int i = 0; i < numHydrophobicBonds; i++) {
      double minValFiveRows = 99999;
      double maxValFiveRows = 0;
      for (int j = 0; j < 5; j++) {//check all 5 remaining motions of the hydrophobic bond
        gsl_matrix_get_row(currentHydrophobicBondRow, hydrophobicBondNullspace, i * 5 + j);
        gsl_vector_minmax(currentHydrophobicBondRow, &minVal, &maxVal);
        if (minVal < minValFiveRows)
          minValFiveRows = minVal;
        if (maxVal > maxValFiveRows)
          maxValFiveRows = maxVal;
      }

      if (minValFiveRows > -RIGID_TOL && maxValFiveRows < RIGID_TOL) {
        gsl_vector_set(m_rigidHydrophobicBonds, i, 1);
        numRigidHydrophobicBonds++;
      }
    }
    log("constraints") << "There are " << numRigidHydrophobicBonds << " rigid out of " << numHydrophobicBonds
                       << " hydrophobic bonds!" << endl;
    gsl_vector_free(currentHydrophobicBondRow);
    gsl_matrix_free(hydrophobicBondNullspace);
  }
}

void Nullspace::projectOnNullSpace(gsl_vector *to_project, gsl_vector *after_project) const {

  if(m_nullspaceSize==0){
    gsl_vector_set_zero(after_project);
    return;
  }

  //The projection equation is:   after_project = N * N^T * to_project
  //Dimensions are: [n x 1] = [n x (n-r) ] [(n-r) x n] [n x 1]
  //Computations: N^T * to_project:	(n-r) * n  --> firstResult [(n-r) x 1]
  //Computations: N * firstResult:	(n) * (n-r)  --> finalResult [n x 1]
  //Computations overall: 2*n*(n-r)

  gsl_vector* firstResult = gsl_vector_alloc(m_nullspaceSize);
  gsl_blas_dgemv (CblasTrans, 1.0, m_nullspaceBasis, to_project, 0.0, firstResult);
  gsl_blas_dgemv (CblasNoTrans, 1.0, m_nullspaceBasis, firstResult, 0.0, after_project);
  gsl_vector_free(firstResult);


  //Old version
//	gsl_blas_dgemv (CblasNoTrans, 1.0, nullspace->P, to_project, 0.0, after_project);
}


gsl_matrix *Nullspace::getBasis() const {
  return m_nullspaceBasis;
}




