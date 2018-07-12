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
#include <iostream>
#include <list>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_vector.h>
#include <math/gsl_helpers.h>
#include <math/NullspaceSVD.h>

#include "core/Molecule.h"
#include "core/Grid.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "applications/options/RigidityOptions.h"
#include "../core/Configuration.h"

extern double jacobianAndNullspaceTime;
extern double rigidityTime;

using namespace std;

int main( int argc, char* argv[] ){

  enableLogger("rigidity");

  if(argc<2){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}

  RigidityOptions::createOptions(argc,argv);
  RigidityOptions& options = *(RigidityOptions::getOptions());

  options.print();

  string out_path = options.workingDirectory;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einfluss von Hydrophobics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ToDo
  // 1. Code-Ablauf der Rigidity Analysis nachvollziehen
  // 2. Quasi-identische Implementierung für Neue Klasse HydrophobicBonds (vgl. mit HBonds, bereits angefangen)
  // 3. Anpassung der Jacobi-Matrix zum Handling mehrerer Constraint-Typen (vgl. mit DBonds, bereits angefangen)
  // 4. Anpassung der Matrix, die freie Bewegungen der Constraint-Typen überprüft (HBondJacobian)
  // 5. Evtl extrahieren der Berechnung der Ableitung in eigene Constraint-Klasse
  // 6. Berechnung des Nullraums
  // 7. Identifikation starrer dihedrals und hydrophobic bonds (vgl. bisher HBonds)

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Einfluss von Hydrophobics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      options.extraCovBonds,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  protein->initializeTree(movingResidues,1.0,options.roots);
  string name = protein->getName();

  log("rigidity")<<"Molecule has:"<<endl;
  log("rigidity") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("rigidity")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("rigidity")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" total bond constraints"<<endl;
  log("rigidity")<<"> "<<protein->getHBonds().size()<<" hydrogen bonds"<<endl;
  log("rigidity")<<"> "<<protein->getHydrophobicBonds().size()<<" hydrophobic bonds"<<endl;
  log("rigidity") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  Configuration* conf = protein->m_conf;
  NullspaceSVD ns = *(dynamic_cast<NullspaceSVD*>(conf->getNullspace()));
  int numRows = ns.getMatrix()->size1;
  int numCols = ns.getMatrix()->size2;
  int nullspaceCols = ns.getNullspaceSize();
  int rankJacobian=min(numCols,numRows)-nullspaceCols;
  int numredundant=numRows-rankJacobian;

  log("rigidity") << "Dimension of Jacobian: " << numRows << " rows, ";
  log("rigidity") << numCols << " columns" << endl;
  log("rigidity") << "Dimension of kernel: " << nullspaceCols << endl;
  log("rigidity") << "Rank of Jacobian: " <<rankJacobian << endl;
  log("rigidity") << "Number of redundant constraints:" <<numredundant << endl;

  /// Create larger rigid substructures for rigid cluster decomposition
  Molecule* rigidified = protein->collapseRigidBonds(2);

  ///Write PDB File for pyMol usage
  int sample_id = 1;
  string out_file = out_path + "output/" + name + "_new_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".pdb";

  rigidified->writeRigidbodyIDToBFactor();
  IO::writePdb(rigidified, out_file);

  if(options.saveData <= 0) return 0;

  ///save pyMol coloring script
  string pyMol = out_path + "output/" + name + "_pyMol_" +
                 std::to_string((long long) sample_id)
                 + ".pml";
  string statFile = out_path + "output/" + name + "_stats_" +
                    std::to_string((long long) sample_id)
                    + ".txt";
  ///Write statistics
  IO::writeStats(protein, statFile, rigidified); //original protein with all bonds etc, rigidified one for cluster info

  ///Write pyMol script
  IO::writePyMolScript(rigidified, out_file, pyMol, protein);

  ///save singular values
  //NullspaceSVD* derived = dynamic_cast<NullspaceSVD*>(conf->getNullspace());
  string outSing = out_path + "output/singVals.txt";
  gsl_vector_outtofile(ns.getSVD()->S, outSing);

  if(options.saveData <= 1) return 0;

  ///save Jacobian and Nullspace to file
  string outJac=out_path + "output/" +  name + "_jac_" +
                std::to_string((long long)sample_id)
                + ".txt";
  string outNull=out_path + "output/" +  name + "_nullSpace_" +
                 std::to_string((long long)sample_id)
                 + ".txt";

  ns.writeMatricesToFiles(outJac, outNull);

  if(options.saveData <= 2) return 0;

  if (conf->getHydrogenJacobian()){
    string outHydrogenJacobian=out_path + "output/" +  name + "_HBondJacobian_" +
                               std::to_string((long long)sample_id)
                               + ".txt";
    gsl_matrix_outtofile(conf->getHydrogenJacobian(), outHydrogenJacobian);
  }

  if (conf->getHydrophobicJacobian()) {
    string outHpJacobian = out_path + "output/" + name + "_HydrophobicBondJacobian_" +
                           std::to_string((long long) sample_id)
                           + ".txt";
    gsl_matrix_outtofile(conf->getHydrophobicJacobian(), outHpJacobian);
  }

  string rbFile=out_path + "output/" +  name + "_RBs_" +
                std::to_string((long long)sample_id)
                + ".txt";
  string covFile=out_path + "output/" +  name + "_covBonds_" +
                 std::to_string((long long)sample_id)
                 + ".txt";

  ///Write covalent bonds
  IO::writeCovBonds(protein,covFile);

  ///Write rigid bodies
  IO::writeRBs(rigidified, rbFile);

  return 0;
}


