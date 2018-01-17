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

extern double jacobianAndNullspaceTime;
extern double rigidityTime;

using namespace std;

int main( int argc, char* argv[] ){

  enableLogger("rigidity");

  if(argc<2){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}

  //RigidityOptions options(argc,argv);
  RigidityOptions::createOptions(argc,argv);
  RigidityOptions& options = *(RigidityOptions::getOptions());

  options.print();

  string out_path = options.workingDirectory;
  //string pdb_file = path + protein_name + ".pdb";

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
  cout<<"Here 1"<<endl;
  protein->initializeTree(movingResidues,1.0,options.roots);
  string name = protein->getName();

  log("rigidity")<<"Molecule has:"<<endl;
  log("rigidity") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("rigidity")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("rigidity")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" bond constraints"<<endl;
  log("rigidity") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  Configuration* conf = protein->m_conf;
  NullspaceSVD ns = *(dynamic_cast<NullspaceSVD*>(conf->getNullspace()));
  int numCols = ns.getMatrix()->size2;
  int nullspaceCols = ns.getNullspaceSize();
  int sampleCount = 0;
  gsl_vector* singValVector = gsl_vector_copy(ns.getSVD()->S);

  log("rigidity") << "Dimension of Jacobian: " << ns.getMatrix()->size1 << " rows, ";
  log("rigidity") << numCols << " columns" << endl;
  log("rigidity") << "Dimension of kernel " << nullspaceCols << endl;

//  if (5*ns.getMatrix()->size1 < numCols){//less constraints than cycle-dofs
//    //Write the complete J*V product out to file
//    gsl_matrix* fullProduct = gsl_matrix_alloc(baseJacobian->size1, numCols);
//    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, baseJacobian, baseNullspaceV, 0.0, fullProduct);
//    if(options.saveData > 1) {
//      string outProd = "fullProduct_JV.txt";
//      gsl_matrix_outtofile(fullProduct, outProd);
//      gsl_matrix_free(fullProduct);
//
//      string outMat = "Vmatrix.txt";
//      gsl_matrix_outtofile(baseNullspaceV, outMat);
//    }
//  }



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
  NullspaceSVD* derived = dynamic_cast<NullspaceSVD*>(conf->getNullspace());
  if (derived) {
    string outSing = out_path + "output/singVals.txt";
    gsl_vector_outtofile(derived->getSVD()->S, outSing);
  }

  if(options.saveData <= 1) return 0;

  ///save Jacobian and Nullspace to file
  string outJac=out_path + "output/" +  name + "_jac_" +
                std::to_string((long long)sample_id)
                + ".txt";
  string outNull=out_path + "output/" +  name + "_nullSpace_" +
                 std::to_string((long long)sample_id)
                 + ".txt";

  if (derived) {
    derived->writeMatricesToFiles(outJac, outNull);
  }

  if(options.saveData <= 2) return 0;

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


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATE ROUCHE-CAPELLI-THEOREM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ToDo
  // 1. Target Protein Struktur einlesen; dazu die "RigidityOptions" erweitern um Target-Struktur
  // 2. gleiches Set an H-bonds (Constraints) für Protein und Target erstellen
  // 3. Konfiguration auf Target erzeugen
  // 4. berechne Rotationswinkel q_init und q_target, berechne Delta_q, "siehe global torsions"
  // 5. Evaluiere Rouche-Capelli Theorem
  // 6. Berechne Störung P einzelner H-bonds mittels S = J * Delta_q
  // 7. Identifiziere Constraints mit großem S
  // 8. Identisches Vorgehen für Delta_p anstelle Delta_q (Wiederholung ab 4.)

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVALUATE ROUCHE-CAPELLI-THEOREM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  return 0;
}


