
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>

#include "Nullspace.h"
#include "Logger.h"
#include "gsl_helpers.h"

using namespace std;

Nullspace::Nullspace(SVD * svd) :
    svd(svd),
    m(svd->matrix->size1),
    n(svd->matrix->size2),
    numCoordinatedDihedrals( 0 ),
    numRigidDihedrals( 0 ),
    numRigidHBonds( 0 ),
    m_nullspaceBasis(nullptr),
    rigidAngles(gsl_vector_alloc(n)),
    rigidHBonds(gsl_vector_alloc(n))
{
}

Nullspace::~Nullspace ()
{
  if(m_nullspaceBasis)
    gsl_matrix_free(m_nullspaceBasis);

  gsl_vector_free(rigidAngles);
  gsl_vector_free(rigidHBonds);
  //delete svd;
}

void Nullspace::UpdateFromMatrix()
{
  svd->UpdateFromMatrix();

  //Compute nullspacesize
  double maxSingularValue = gsl_vector_get(svd->S,0);

  //Case with m < n and all singular values non-zero
  m_nullspaceSize = std::max( (int)(svd->V->size2 - svd->matrix->size1),0);

  for (int i=0; i<svd->S->size ; ++i){
    if ( gsl_vector_get(svd->S,i) / maxSingularValue < SINGVAL_TOL ) {
      m_nullspaceSize = svd->V->size2-i;
      break;
    }
  }

  //TODO: If an existing basis of proper size is allocated we might not need to reallocate here
  if(m_nullspaceBasis)
    gsl_matrix_free(m_nullspaceBasis);

  if(m_nullspaceSize > 0) {
    gsl_matrix_view nullspaceBasis_view = gsl_matrix_submatrix(svd->V,
                                                               0,                               //Row
                                                               svd->V->size2 - m_nullspaceSize, //Col
                                                               svd->V->size2,                 //Height
                                                               m_nullspaceSize);                  //Width

    m_nullspaceBasis = gsl_matrix_calloc(svd->V->size2, m_nullspaceSize);
    gsl_matrix_memcpy(m_nullspaceBasis, &nullspaceBasis_view.matrix);
  }
  else{
    m_nullspaceBasis = gsl_matrix_calloc(svd->V->size2,1);//1-dim vector with zeros as entries
  }
}


/** Analyzes which dihedrals and hydrogen bonds are rigidified by constraints */
void Nullspace::RigidityAnalysis(gsl_matrix* HBondJacobian)
{
  // First, check the dihedral angles for rigidity

  // Nullspace dimension
  gsl_vector* currentRow = gsl_vector_alloc(n);

  gsl_vector_set_zero(rigidAngles);
  gsl_vector_set_zero(rigidHBonds);

  bool moving=false;
  numCoordinatedDihedrals = 0;
  numRigidDihedrals = 0;

  for(int i=0; i<n; i++){
    gsl_matrix_get_row(currentRow,svd->V,i);

    for(int j=n-m_nullspaceSize; j<n; j++){
      double val = fabs( gsl_vector_get(currentRow,j) );
      if( val > RIGID_TOL ){
        moving = true;
      }
    }

    // at least one entry was greater than threshold --> dihedral is coordinated
    if( moving ) {
      numCoordinatedDihedrals++;
    }
    else {
      gsl_vector_set(rigidAngles,i,1); /// binary list of rigid dihedrals, complement to coordinated version
      numRigidDihedrals++;
    }
    // reset for the next dihedral / row
    moving = false;
  }

  log("constraints")<<"There are "<<numRigidDihedrals << " rigidified and " << numCoordinatedDihedrals << " coordinated dihedrals" << endl;

  // Free memory
  gsl_vector_free(currentRow);

  // Now, check the hydrogen Bonds for rigidity
  int numHBonds = HBondJacobian->size1;
  gsl_matrix* hBondNullspace;// = gsl_matrix_alloc(numHBonds,nullspaceSize);
  gsl_vector* currentHBondRow;// = gsl_vector_alloc(nullspaceSize);
  if(m_nullspaceSize==0){
    hBondNullspace = gsl_matrix_calloc(numHBonds,1);
    currentHBondRow = gsl_vector_calloc(1);
  }else{
    hBondNullspace = gsl_matrix_alloc(numHBonds,m_nullspaceSize);
    currentHBondRow = gsl_vector_alloc(m_nullspaceSize);
  }

  ///Calculate the "In-Nullspace" Rotation of the hBonds
  //gsl_matrix_view m_nullspaceBasis = gsl_matrix_submatrix(svd->V,0,svd->V->size2-m_nullspaceSize,n,m_nullspaceSize); //Matrix, top-left element, num rows, num cols
  //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HBondJacobian, &m_nullspaceBasis.matrix, 0.0, hBondNullspace);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, HBondJacobian, m_nullspaceBasis, 0.0, hBondNullspace);

  //Only for test analysis: write hBondNullspace=J_h*N to file
//	string outNullH="KINARI_Comparison/output/hBondnullSpace.txt";
//	gsl_matrix_outtofile(hBondNullspace, outNullH);

  double maxVal;
  double minVal;
  numRigidHBonds = 0;

  for(int i=0; i<numHBonds; i++){
    gsl_matrix_get_row(currentHBondRow,hBondNullspace,i);
    gsl_vector_minmax(currentHBondRow, &minVal, &maxVal);

    if (minVal > -RIGID_TOL && maxVal < RIGID_TOL){
      gsl_vector_set(rigidHBonds,i,1);
      numRigidHBonds++;
    }
  }

//	log("constraints")<<"There are "<<m_numRigidHBonds<<" rigid out of "<<numHBonds<<" hydrogen bonds!"<<endl;
	log("constraints")<<"There are "<<numRigidHBonds<<" rigid out of "<<numHBonds<<" hydrogen bonds!"<<endl;

  gsl_vector_free(currentHBondRow);
  gsl_matrix_free(hBondNullspace);
}

void Nullspace::ProjectOnNullSpace (gsl_vector *to_project, gsl_vector *after_project) const {

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


void Nullspace::WriteMatricesToFiles(
    const std::string& jac_file,
    const std::string& null_file,
    const std::string& sval_file) const
{
  gsl_matrix_outtofile(svd->matrix, jac_file);
  gsl_matrix_outtofile(svd->V, null_file);
  gsl_vector_outtofile(svd->S, sval_file);
}

gsl_matrix *Nullspace::getBasis() const {
  return m_nullspaceBasis;
}

SVD *Nullspace::getSVD() const {
  return svd;
}

