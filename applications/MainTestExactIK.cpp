//
// Created by Rasmus Fonseca on 29/03/16.
//

#include <string>
#include <Logger.h>
#include <math/MKLSVD.h>
#include <math/GSLSVD.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "core/Chain.h"
#include "IO.h"
#include "loopclosure/ExactIK.h"
#include "SamplingOptions.h"
#include "math/gsl_helpers.h"

int main( int argc, char* argv[] ) {

  string pdb_file = "/Users/rfonseca/Downloads/1crn_Processed.pdb";
  Molecule protein;
  vector<string> extracovbonds;
  IO::readPdb( &protein, pdb_file, extracovbonds );

  IO::readHbonds( &protein, "/Users/rfonseca/Downloads/hbonds.txt" );

  //Create the rigid body trees
  IO::readRigidbody( &protein );

  protein.buildRigidbodyTree(0, false);//with the rigid body tree in place, we can generate a configuration

  ExactIK ik;

  int first = 30;
  vector<Configuration*> rebuiltConfs = ik.rebuildLoop(
      protein.getChain("A")->getResidue(first+0),
      protein.getChain("A")->getResidue(first+1),
      protein.getChain("A")->getResidue(first+5)
  );

  cout<<rebuiltConfs.size()<<" rebuilt configurations"<<endl;
  for(size_t i=0;i<rebuiltConfs.size(); i++){
    IO::writePdb(rebuiltConfs[i]->updatedProtein(), "/Users/rfonseca/Downloads/rebuilt_"+std::to_string(i)+".pdb");
  }
}
