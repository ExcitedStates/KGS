#include <string>
#include <Logger.h>
#include <math/MKLSVD.h>
#include <math/GSLSVD.h>
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

int main( int argc, char* argv[] ) {
  std::string prg = argv[0];
  argc=9;
  argv = new char *[argc];
  argv[0] = "test";
  argv[1] = "--initial";
  argv[2] = "/Users/rfonseca/Downloads/1crn_Processed.pdb";
  argv[3] = "--hbondMethod";
  argv[4] = "user";
  argv[5] = "--hbondFile";
  argv[6] = "/Users/rfonseca/Downloads/hbonds.txt";
  argv[7] = "--collisionFactor";
  argv[8] = "0.0";

  SamplingOptions::createOptions(argc, argv);
  Molecule protein;
  vector<string> extracovbonds;
  IO::readPdb( &protein, SamplingOptions::getOptions()->initialStructureFile, extracovbonds );

  IO::readHbonds( &protein, "/Users/rfonseca/Downloads/hbonds.txt" );

  //Create the rigid body trees
  IO::readRigidbody( &protein );
  protein.buildSpanningTree(0, false);

  Configuration* conf = new Configuration(&protein);

  IO::writePdb(conf->updatedProtein(), "/Users/rfonseca/Downloads/test1.pdb");

  gsl_vector* gradient = gsl_vector_calloc(protein.m_spanning_tree->m_numDOFs);
  gsl_vector_set(gradient, 4, 0.1);
  Move* move = new RawMove();
  for(int i=0;i<10;i++){
    Configuration* cNew = move->move(conf, gradient);
    if(cNew==NULL){
      cerr<<"ERROR .. new conf is null"<<endl;
      continue;
    }
    IO::writePdb(cNew->updatedProtein(), "/Users/rfonseca/Downloads/test"+std::to_string(i)+".pdb");
    conf = cNew;

  }

  std::system("/usr/local/bin/pymol /Users/rfonseca/Downloads/test*pdb");


}
