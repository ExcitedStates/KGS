

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
#include "CTKTimer.h"

using namespace std;

extern double prevRMSDTime;
extern double currentRMSDTime;

Molecule * myReadFile(string pdbFile){
  char* tmp = realpath(pdbFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<pdbFile<<" is not a valid PDB-file"<<endl; exit(-1); }

  Selection movingResidues("all");
  Molecule* protein = IO::readPdb(
      pdbFile,
      movingResidues
  );
  protein->setCollisionFactor(1.0);

  return protein;
}


int main( int argc, char* argv[] ){
  enableLogger("rmsd");
  if(argc<3){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}
  
  Selection sel("all");
  Configuration* reference = new Configuration(myReadFile(argv[1]));
  for(int i=2;i<argc;i++){
    Molecule * p = myReadFile(argv[i]);
    Configuration* c = new Configuration(p);

    CTKTimer timer;
    timer.Reset();
    double start_time = timer.LastElapsedTime();

    double rmsd = p->alignReferencePositionsTo(reference->getMolecule(), sel);
//        metric->align(p, reference->getMolecule());
    log("rmsd")<<argv[i]<<" : "<<rmsd<<endl;

    double end_time = timer.ElapsedTime();

    string outfile = p->getName()+"_aligned.pdb";
    IO::writePdb(p,outfile);

    cout<<"Alignment took "<<end_time-start_time<<endl;
//    cout<<"RMSD only took "<<end_time_2-end_time<<endl;
    delete c;
    delete p;
  }

}


