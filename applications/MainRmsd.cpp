

#include <string>
#include <iostream>
#include <metrics/Dihedral.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "metrics/RMSD.h"

using namespace std;

Molecule * myReadFile(string pdbFile){
    char* tmp = realpath(pdbFile.c_str(), nullptr);
    if(tmp==nullptr){ cerr<<pdbFile<<" is not a valid PDB-file"<<endl; exit(-1); }
    Molecule * protein = new Molecule();
    vector<string> extraCovBonds;
    IO::readPdb( protein, pdbFile, extraCovBonds );
    IO::readRigidbody( protein );
    protein->buildSpanningTree();

    return protein;
}


int main( int argc, char* argv[] ){
  enableLogger("rmsd");
  if(argc<3){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}

//  metrics::Metric* metric = new metrics::RMSD();
  metrics::Metric* metric = new metrics::Dihedral();
  Configuration* reference = new Configuration(myReadFile(argv[1]));
  for(int i=2;i<argc;i++){
    Molecule * p = myReadFile(argv[i]);
    Configuration* c = new Configuration(p);
    double rmsd = metric->distance(c, reference);
    log("rmsd")<<argv[i]<<" : "<<rmsd<<endl;
    delete p;
    delete c;
  }

}


