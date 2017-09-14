/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

This entire text, including the above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

*/


#include <string>
#include <vector>
#include <ctime>
#include <Logger.h>
#include <math/SVDMKL.h>
#include <math/SVDGSL.h>
#include <math/CudaSVD.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <moves/Move.h>
#include <moves/NullspaceMove.h>
#include <moves/RawMove.h>
#include <directions/Direction.h>
#include <directions/MSDDirection.h>
#include <metrics/RMSDnosuper.h>
#include <directions/BlendedDirection.h>
#include <directions/RandomDirection.h>
#include <sys/time.h>
#include <math/QR.h>
#include <math/QRMKL.h>
#include <math/QRGSL.h>
#include <math/NullspaceSVD.h>
#include <math/NullspaceQR.h>

#include "core/Chain.h"
#include "IO.h"
#include "loopclosure/ExactIK.h"
#include "applications/options/ExploreOptions.h"
#include "math/gsl_helpers.h"
#include "core/Configuration.h"

using namespace std;

void testGlobalMSD();
void testGlobalGradient();
void testQR();
void testSelection();
void testIncDOFs();

int main( int argc, char* argv[] ) {
  enableLogger("default");
//  testGlobalMSD();
//  testGlobalGradient();
//  testQR();
//  testSelection();
  testIncDOFs();
}


void testIncDOFs(){
  Selection sel("all");
  Molecule* mol = IO::readPdb("/Users/rfonseca/helix.pdb", sel );
  int n = mol->m_spanningTree->getNumDOFs();

  Configuration* finalconf = new Configuration(mol);
  for(int d=6;d<n;d++){
    finalconf->m_dofs[d] = 50*3.141592/180.0;
  }
  IO::writePdb(finalconf->updatedMolecule(), "/Users/rfonseca/helix_out_final.pdb");

  mol->setConfiguration(new Configuration(mol));

  for(int i=0;i<50;i++){
    Configuration* conf = new Configuration(mol->m_conf);
    for(int d=6;d<n;d++){
      conf->m_dofs[d] = 1*3.141592/180.0;
    }
    IO::writePdb(conf->updatedMolecule(), "/Users/rfonseca/helix_out_"+to_string(i)+".pdb");

    for(Atom* a: mol->getAtoms()){
      a->m_referencePosition = a->m_position;
    }
  }


}

//From http://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

using namespace metrics;

void testSelection(){
  try {
    enableLogger("debug");
    Selection moving("all");
    std::vector<std::string> extraCovBonds;
    Molecule* mol = IO::readPdb("/Users/rfonseca/2kmj_1.pdb", moving, extraCovBonds);
    Selection sel("resi 17 and id 4+7");
    for(auto a: sel.getSelectedAtoms(mol)){
      std::cout<<a<<endl;
    }

  } catch (const std::string& ex) {
    cerr<<ex<<endl;
  }
}
void testGlobalGradient(){
  try {
    enableLogger("debug");
    std::vector<std::string> extraCovBonds;
    Selection sel("all");
    Molecule* mol = IO::readPdb("/Users/rfonseca/Y.pdb", sel, extraCovBonds);

    Configuration* conf = new Configuration(mol);
    cout<<"DOFS: "<<conf->getNumDOFs()<<endl;

    int dof = 8;
    Coordinate pos = conf->updatedMolecule()->getAtom("A", 29, "OH")->m_position;
    cout<<"Old OH pos: "<<pos<<endl;
    Math3D::Vector3 der = mol->m_spanningTree->Edges.at(dof)->getDOF()->getDerivative(pos);
    der = der*0.01;
    cout<<der<<endl;
    Configuration* conf2 = new Configuration(conf);
    conf2->m_dofs[dof]+=0.01;
    Coordinate pos2 = conf2->updatedMolecule()->getAtom("A", 29, "OH")->m_position;
    cout<<"New OH pos: "<<pos2<<endl;
    cout<<"Expected:   "<<(pos+der)<<endl;
    if(mol->m_spanningTree->Edges.at(dof)->getBond()==nullptr) {
      cout << "Global DOF" << endl;
    }else {
      cout << "DOF bond:   " << (*mol->m_spanningTree->Edges.at(dof)->getBond()) << endl;
    }

  } catch (const std::string& ex) {
    cerr<<ex<<endl;
  }
}


void testGetTorsion(){
  try {
    Selection all("all");
    std::vector<std::string> extraCovBonds;
    Molecule* mol = IO::readPdb("/Users/rfonseca/1crn.pdb", all, extraCovBonds);
    Selection sel("backbone");
    for (auto const &b: sel.getSelectedBonds(mol)) {
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
    for(int j=0;j<n;j++){
      if(rand()%100<10)
        gsl_matrix_set(A, i, j, (rand()*4.0/RAND_MAX));
    }

//  gsl_matrix_cout(A);

  double start = get_wall_time();
  SVD* mklsvd = new SVDMKL(A);
  mklsvd->UpdateFromMatrix();
//  mklsvd->print();
  double duration = ( get_wall_time() - start );
  cout<<"MKL took  "<<duration<<"secs"<<endl;


//  gsl_matrix* S = gsl_matrix_calloc(A->size1, A->size2);
//  for(int i=0;i<std::min(A->size1,A->size2);i++){
//    gsl_matrix_set(S, i,i, gsl_vector_get(mklsvd->S, i));
//  }
//  cout<<"Product: "<<endl;
//  gsl_matrix_cout( gsl_matrix_mul(mklsvd->U, gsl_matrix_mul(S, gsl_matrix_trans(mklsvd->V))) );

//  cout<<" ------------------------------------- "<<endl;

  start = get_wall_time();
  SVD* cudasvd = new CudaSVD(A);
  cudasvd->UpdateFromMatrix();
//  cudasvd->print();
  duration = ( get_wall_time() - start );
  cout<<"CUDA took "<<duration<<"secs"<<endl;

  start = get_wall_time();
  SVD* gslsvd = new SVDGSL(A);
  gslsvd->UpdateFromMatrix();
  duration = ( get_wall_time() - start );
  cout<<"GSL took  "<<duration<<"secs"<<endl;

//  S = gsl_matrix_calloc(A->size1, A->size2);
//  for(int i=0;i<std::min(A->size1,A->size2);i++){
//    gsl_matrix_set(S, i,i, gsl_vector_get(cudasvd->S, i));
//  }
//
//  cout<<"Product: "<<endl;
//  gsl_matrix_cout( gsl_matrix_mul(cudasvd->U, gsl_matrix_mul(S, gsl_matrix_trans(cudasvd->V))) );

}

//void testQR(){
////  srand(101);
////  int m = 4;
////  int n = 3;
////  gsl_matrix* A  = gsl_matrix_calloc(m,n);
////  rand();rand();rand();
////  for(int i=0;i<m;i++) {
////    for (int j = 0; j < n; j++) {
////      gsl_matrix_set(A, i, j, rand() % 4);
////    }
////  }
//  srand(14);
//  rand();rand();rand();
////  int m = 7;
////  int n = 10;
//  int m = 5000;
//  int n = 1000;
//  cout<<m<<" x "<<n<<endl;
//  gsl_matrix* A  = gsl_matrix_calloc(m,n);
//  for(int i=0;i<m;i++) {
//    for (int j = 0; j < n; j++) {
//      gsl_matrix_set(A, i, j, (rand()%2==0)?(rand()*1.0/RAND_MAX):(0.0));
////      gsl_matrix_set(A, i, j, rand()%3);
//    }
//  }
//
//  for(int j=0;j<n;j++){
//    gsl_matrix_set(A,3,j, gsl_matrix_get(A,2,j));
//  }
//
////  std::cout<<"A:"<<endl;
////  for(int i=0;i<m;i++){ for(int j=0;j<n;j++) printf(" %5.2f",gsl_matrix_get(A,i,j)); std::cout<<std::endl;}
//
//
//  double start, duration;
//
//  start = get_wall_time();
//  SVD* svd = new SVDMKL(A);//SVD::createSVD(A);
//  Nullspace* ns2 = new NullspaceSVD(svd);
//  ns2->updateFromMatrix();
//  duration = ( get_wall_time() - start );
//  cout<<"SVD took "<<duration<<"secs"<<endl;
////  std::cout<<"N_SVD:"<<endl;
////  gsl_matrix_cout(ns2.getBasis());
//
//  start = get_wall_time();
//  TransposeQR* qr = new TransposeQR(A);
//  Nullspace* ns = new NullspaceQR(qr);
//  ns->updateFromMatrix();
//  duration = ( get_wall_time() - start );
//  cout<<"QR took "<<duration<<"secs"<<endl;
////  std::cout<<"N:"<<endl;
////  gsl_matrix_cout(ns.getBasis());
//
////  std::cout<<"A·N:"<<endl;
////  gsl_matrix_cout(gsl_matrix_mul(A, ns.getBasis()));
////  std::cout<<"A·N_SVD:"<<endl;
////  gsl_matrix_cout(gsl_matrix_mul(A, ns2.getBasis()));
//
//
////  gsl_matrix_cout(A);
////
////  double start, duration;
//
////  start = get_wall_time();
////  QR* mklqr = new QRMKL(A);
////  mklqr->updateFromMatrix();
////  mklqr->print();
////  duration = ( get_wall_time() - start );
////  cout<<"MKL took  "<<duration<<"secs"<<endl;
//
//
////  gsl_matrix* S = gsl_matrix_calloc(A->size1, A->size2);
////  for(int i=0;i<std::min(A->size1,A->size2);i++){
////    gsl_matrix_set(S, i,i, gsl_vector_get(mklqr->S, i));
////  }
////  cout<<"Product: "<<endl;
////  gsl_matrix_cout( gsl_matrix_mul(mklqr->U, gsl_matrix_mul(S, gsl_matrix_trans(mklqr->V))) );
//
////  cout<<" ------------------------------------- "<<endl;
//
////  start = get_wall_time();
////  QR* gslqr = new QRGSL(A);
////  gslqr->updateFromMatrix();
////  gslqr->print();
////  duration = ( get_wall_time() - start );
////  cout<<"GSL took  "<<duration<<"secs"<<endl;
//
////  S = gsl_matrix_calloc(A->size1, A->size2);
////  for(int i=0;i<std::min(A->size1,A->size2);i++){
////    gsl_matrix_set(S, i,i, gsl_vector_get(cudasvd->S, i));
////  }
////
////  cout<<"Product: "<<endl;
////  gsl_matrix_cout( gsl_matrix_mul(cudasvd->U, gsl_matrix_mul(S, gsl_matrix_trans(cudasvd->V))) );
//
//}
