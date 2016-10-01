//
// Created by Rasmus Fonseca on 29/03/16.
//

#include <string>
#include <Logger.h>
#include <math/SVDMKL.h>
#include <math/SVDGSL.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include "core/Chain.h"
#include "IO.h"
#include "loopclosure/ExactIK.h"
#include "SamplingOptions.h"
#include "math/gsl_helpers.h"

using namespace std;

int main( int argc, char* argv[] ) {

  string pdb_file = "/Users/rfonseca/Documents/KGSexperiments/WAFR16/work/dhfr/MET20+FG_rigidCore_ligands/3ql3_clean.Processed.pdb.knr";
  Molecule protein;
  vector<string> extra_cov_bonds;
  IO::readPdb(&protein, pdb_file, extra_cov_bonds);
  protein.setCollisionFactor(0.4);

  IO::readHbonds(&protein,
                 "/Users/rfonseca/Documents/KGSexperiments/WAFR16/work/dhfr/MET20+FG_rigidCore_ligands/constraints.txt");

  //Create the rigid body trees
  IO::readRigidbody(&protein);

  protein.buildSpanningTree();//with the rigid body tree in place, we can generate a configuration
  Configuration* init = new Configuration(&protein);
  protein.setConfiguration(init);
  protein.m_initialCollisions = protein.getAllCollisions();

  ExactIK ik;

  vector<pair<int, int> > intervals = {{20,  25},
                                       {116, 121},
                                       {124, 128}};

  vector< tuple<int,int,int> > triples;

  //Generate all triples within intervals
  for( auto ival: intervals ) {
    for( int r1=ival.first; r1<=ival.second; r1++ ) {
      for( int r2=r1 + 1; r2<=ival.second; r2++ ) {
        for( int r3=r2 + 1; r3<=ival.second; r3++ ) {
          triples.push_back(make_tuple(r1, r2, r3));
        }
      }
    }
  }


  vector<Configuration *> allConfs;

  for (auto ival: intervals) {
    for (int r1 = ival.first; r1 <= ival.second; r1++) {
      for (int r2 = r1 + 1; r2 <= ival.second; r2++) {
        for (int r3 = r2 + 1; r3 <= ival.second; r3++) {

          Residue* res1 = protein.getChain("A")->getResidue(r1);
          Residue* res2 = protein.getChain("A")->getResidue(r2);
          Residue* res3 = protein.getChain("A")->getResidue(r3);
          if(!ik.validRebuildLoop(res1,res2,res3)) continue;

          protein.setConfiguration(init);
          vector<Configuration *> rebuiltConfs = ik.rebuildLoop(res1,res2,res3);
          int count=0;
          for(auto conf: rebuiltConfs) {
            if(!conf->updatedMolecule()->inCollision()) {
//            if(true) {
              allConfs.push_back(conf);
              count++;
            }

          }
          cout<<"Generated "<<rebuiltConfs.size()<<" configurations ("<<count<<" non-colliding) for triple {"<<r1<<","<<r2<<","<<r3<<"}"<<endl;
        }
      }
    }
  }

  cout<<"Generated "<<allConfs.size()<<" in total"<<endl;
  for (size_t i = 0; i < allConfs.size(); i++) {
    IO::writePdb(allConfs[i]->updatedMolecule(),
                 "/Users/rfonseca/Documents/KGSexperiments/WAFR16/work/dhfr/MET20+FG_rigidCore_ligands_exactik/seeds/rebuilt_" + std::to_string((long long) i) + ".pdb");
  }

//  int first = 30;
//  vector<Configuration *> rebuiltConfs = ik.rebuildLoop(
//      protein.getChain("A")->getResidue(first + 0),
//      protein.getChain("A")->getResidue(first + 1),
//      protein.getChain("A")->getResidue(first + 5)
//  );
//
//  cout << rebuiltConfs.size() << " rebuilt configurations" << endl;
//  for (size_t i = 0; i < rebuiltConfs.size(); i++) {
//    IO::writePdb(rebuiltConfs[i]->updatedMolecule(),
//                 "/Users/rfonseca/Downloads/rebuilt_" + std::to_string((long long) i) + ".pdb");
//  }
}
