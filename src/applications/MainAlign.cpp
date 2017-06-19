

#include <string>
#include <iostream>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "IO.h"
#include "Logger.h"
#include "metrics/RMSD.h"
#include "CTKTimer.h"
#include "applications/options/AlignOptions.h"

using namespace std;

extern double prevRMSDTime;
extern double currentRMSDTime;

Molecule * myReadFile(string pdbFile){
  char* tmp = realpath(pdbFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<pdbFile<<" is not a valid PDB-file"<<endl; exit(-1); }

  Molecule* protein = IO::readPdb(pdbFile);
  protein->setCollisionFactor(1.0);

  return protein;
}


int main( int argc, char* argv[] ){
  enableLogger("align");

  AlignOptions::createOptions(argc, argv);
  AlignOptions &options = *(AlignOptions::getOptions());

  Selection alignSel(options.alignSelection);
  Configuration* reference = new Configuration(myReadFile(options.initialStructureFile));

  Molecule * p = myReadFile(options.targetStructureFile);
  Configuration* c = new Configuration(p);

  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  double rmsd = p->alignReferencePositionsTo(reference->getMolecule(), alignSel);
//        metric->align(p, reference->getMolecule());
  log("align")<<options.targetStructureFile<<" : "<<rmsd<<endl;

  double end_time = timer.ElapsedTime();

  string outfile = p->getName()+"_aligned.pdb";
  IO::writePdb(p,outfile);

  cout<<"Alignment took "<<end_time-start_time<<endl;
  delete c;
  delete p;
  delete reference;
}


