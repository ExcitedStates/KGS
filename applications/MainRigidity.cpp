#include <string>
#include <iostream>
#include <list>

#include <gsl/gsl_matrix.h>
#include <math/NullspaceSVD.h>

#include "core/Molecule.h"
#include "core/Grid.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "applications/options/SamplingOptions.h"

extern double jacobianTime;
extern double rigidityTime;

using namespace std;

int main( int argc, char* argv[] ){

  enableLogger("rigidity");

  if(argc<2){ cerr<<"Too few arguments. Please specify PDB-file in arguments"<<endl; exit(-1);}

  //SamplingOptions options(argc,argv);
  SamplingOptions::createOptions(argc,argv);
  SamplingOptions& options = *(SamplingOptions::getOptions());

  string out_path = options.workingDirectory;
  //string pdb_file = path + protein_name + ".pdb";

  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      movingResidues,
      options.extraCovBonds,
      options.roots,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  string name = protein->getName();

//  string hBondIn = "hBonds_in.txt";
//  IO::writeHbondsIn(protein,hBondIn );

  log("rigidity")<<"Molecule has:"<<endl;
  log("rigidity") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("rigidity")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("rigidity")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("rigidity") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  Configuration* conf = protein->m_conf;
//  Configuration* conf = new Configuration(protein);
//  protein->setConfiguration(conf);
  //conf->computeCycleJacobianAndNullSpace();

  log("rigidity")<<"Dimension of Jacobian: " << conf->getNullspace()->getMatrix()->size1 << " rows, ";
  log("rigidity")<< conf->getNullspace()->getMatrix()->size2<<" columns"<<endl;
  log("rigidity")<<"Dimension of kernel "<< conf->getNullspace()->getNullspaceSize()<<endl;

  int sample_id = 1;
  string out_file = out_path + "output/" + name + "_new_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".pdb";

  Molecule* rigidified = protein->collapseRigidBonds(2);

  ///Write PDB File for pyMol usage
  IO::writeRigidbodyIDToBFactor(rigidified);
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
  IO::writeStats(protein, statFile);


  ///Write pyMol script
  IO::writePyMolScript(rigidified, out_file, pyMol);

  if(options.saveData <= 1) return 0;

  ///save Jacobian and Nullspace to file
  string outJac=out_path + "output/" +  name + "_jac_" +
                std::to_string((long long)sample_id)
                + ".txt";
  string outNull=out_path + "output/" +  name + "_nullSpace_" +
                 std::to_string((long long)sample_id)
                 + ".txt";
  ///save singular values
  string outSing=out_path + "output/" +  name + "_singVals_" +
                 std::to_string((long long)sample_id)
                 + ".txt";
  string rbFile=out_path + "output/" +  name + "_RBs_" +
                std::to_string((long long)sample_id)
                + ".txt";
  string covFile=out_path + "output/" +  name + "_covBonds_" +
                 std::to_string((long long)sample_id)
                 + ".txt";

  ///Write Jacobian
  //gsl_matrix_outtofile(m_molecule->m_conf->CycleJacobian, outJac);
  /////Write Null Space matrix
  //gsl_matrix_outtofile(conf->CycleNullSpace->m_nullspaceBasis,outNull);
  /////Write Singular values
  //gsl_vector_outtofile(NullSpaceRet::singularValues,outSing);
  if (NullspaceSVD* derived = dynamic_cast<NullspaceSVD*>(conf->getNullspace())) {
    derived->writeMatricesToFiles(outJac, outNull, outSing);
  }

  ///Write covalent bonds
  IO::writeCovBonds(protein,covFile);

  ///Write rigid bodies
  IO::writeRBs(rigidified, rbFile);

  return 0;
}


