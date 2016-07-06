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

using namespace std;

int main( int argc, char* argv[] ) {

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
