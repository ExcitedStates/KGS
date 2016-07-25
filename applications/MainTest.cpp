#include <string>
#include <ctime>
#include <Logger.h>
#include <math/MKLSVD.h>
#include <math/GSLSVD.h>
#include <math/CudaSVD.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <moves/Move.h>
#include <moves/NullspaceMove.h>
#include <moves/RawMove.h>

#include "core/Chain.h"
#include "IO.h"
#include "loopclosure/ExactIK.h"
#include "SamplingOptions.h"
#include "math/gsl_helpers.h"
#include "core/Configuration.h"

using namespace std;

int main( int argc, char* argv[] ) {
  enableLogger("default");

  srand(101);
  int m = 1000;
  int n = 1000;
  gsl_matrix* A  = gsl_matrix_calloc(m,n);
  rand();rand();rand();
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      gsl_matrix_set(A, i, j, (rand()*4.0/RAND_MAX));

//  gsl_matrix_cout(A);

  time_t start = clock();
  SVD* mklsvd = new MKLSVD(A);
  mklsvd->UpdateFromMatrix();
//  mklsvd->print();
  double duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"MKL took  "<<duration<<"secs"<<endl;


//  gsl_matrix* S = gsl_matrix_calloc(A->size1, A->size2);
//  for(int i=0;i<std::min(A->size1,A->size2);i++){
//    gsl_matrix_set(S, i,i, gsl_vector_get(mklsvd->S, i));
//  }
//  cout<<"Product: "<<endl;
//  gsl_matrix_cout( gsl_matrix_mul(mklsvd->U, gsl_matrix_mul(S, gsl_matrix_trans(mklsvd->V))) );

  cout<<" ------------------------------------- "<<endl;

  start = clock();
  SVD* cudasvd = new CudaSVD(A);
  cudasvd->UpdateFromMatrix();
//  cudasvd->print();
  duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"CUDA took "<<duration<<"secs"<<endl;

//  S = gsl_matrix_calloc(A->size1, A->size2);
//  for(int i=0;i<std::min(A->size1,A->size2);i++){
//    gsl_matrix_set(S, i,i, gsl_vector_get(cudasvd->S, i));
//  }
//
//  cout<<"Product: "<<endl;
//  gsl_matrix_cout( gsl_matrix_mul(cudasvd->U, gsl_matrix_mul(S, gsl_matrix_trans(cudasvd->V))) );
  /*
  try {
    Molecule mol;
    std::vector<std::string> extraCovBonds;
    IO::readPdb(&mol, "/Users/rfonseca/1crn.pdb", extraCovBonds);
    IO::readRigidbody( &mol );
    mol.buildSpanningTree();
    Selection sel("backbone");
    for (auto const &b: sel.getSelectedBonds(&mol)) {
      cout << b << " " << b->getTorsion() <<endl;
    }
  } catch (const std::string& ex) {
    cerr<<ex<<endl;
  }
   */

}
