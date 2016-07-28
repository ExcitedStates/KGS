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

void testGlobalGradient();

int main( int argc, char* argv[] ) {
  enableLogger("default");
  testGlobalGradient();
}


void testGlobalGradient(){
  try {
    enableLogger("debug");
    Molecule* mol = new Molecule();
    std::vector<std::string> extraCovBonds;
    IO::readPdb(mol, "/Users/rfonseca/Y.pdb", extraCovBonds);
    Selection sel("all");
    IO::readRigidbody( mol, sel );
    mol->buildSpanningTree();

    Configuration* conf = new Configuration(mol);
    cout<<"DOFS: "<<conf->getNumDOFs()<<endl;

    int dof = 8;
    Coordinate pos = conf->updatedMolecule()->getAtom("A", 29, "OH")->m_position;
    cout<<"Old OH pos: "<<pos<<endl;
    Math3D::Vector3 der = mol->m_spanning_tree->Edges.at(dof)->getDOF()->getDerivative(pos);
    der = der*0.01;
    cout<<der<<endl;
    Configuration* conf2 = new Configuration(conf);
    conf2->m_dofs[dof]+=0.01;
    Coordinate pos2 = conf2->updatedMolecule()->getAtom("A", 29, "OH")->m_position;
    cout<<"New OH pos: "<<pos2<<endl;
    cout<<"Expected:   "<<(pos+der)<<endl;
    if(mol->m_spanning_tree->Edges.at(dof)->getBond()==nullptr) {
      cout << "Global DOF" << endl;
    }else {
      cout << "DOF bond:   " << (*mol->m_spanning_tree->Edges.at(dof)->getBond()) << endl;
    }

  } catch (const std::string& ex) {
    cerr<<ex<<endl;
  }


}

void testGetTorsion(){
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
}

void testCuda(){
  srand(101);
  int m = 1000;
  int n = 1000;
  gsl_matrix* A  = gsl_matrix_calloc(m,n);
  rand();rand();rand();
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++) {
      if((rand()*1.0/RAND_MAX)<0.1)
        gsl_matrix_set(A, i, j, (rand() * 4.0 / RAND_MAX));
    }

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

}