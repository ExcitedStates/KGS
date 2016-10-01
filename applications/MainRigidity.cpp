#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#include <gsl/gsl_matrix.h>
#include <math/gsl_helpers.h>
#include <gsl/gsl_matrix_double.h>
#include <math/NullspaceSVD.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "planners/SamplingPlanner.h"
#include "planners/RRTPlanner.h"
#include "planners/DihedralRRT.h"
#include "core/Grid.h"
#include "Util.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "core/ProteinHBond.h"
#include "Logger.h"
#include "SamplingOptions.h"
#include "moves/CompositeMove.h"

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

  Molecule * protein = new Molecule();
  IO::readPdb( protein, options.initialStructureFile, options.extraCovBonds );
  string name = protein->getName();

  if(options.hydrogenbondMethod=="user")
    IO::readHbonds( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="rnaview")
    IO::readHbonds_rnaview( protein, options.hydrogenbondFile, options.annotationFile.empty() );
  else if(options.hydrogenbondMethod=="first" || options.hydrogenbondMethod=="FIRST")
    IO::readHbonds_first( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="kinari" || options.hydrogenbondMethod=="KINARI")
    IO::readHbonds_kinari( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="hbplus" || options.hydrogenbondMethod=="hbPlus")
    IO::readHbonds_hbPlus( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="vadar")
    IO::readHbonds_vadar( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="dssr")
    IO::readHbonds_dssr( protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="identify")
    HbondIdentifier::identifyHbonds(protein);

  cout<<"Rigidity: "<<options.hydrogenbondMethod<<endl;

  string hBondIn = "../hBonds_in.txt";
  IO::writeHbondsIn(protein,hBondIn );

  IO::readRigidbody( protein );
  protein->buildSpanningTree();
//	if(options.hydrogenbondMethod!="user")
//		writeHBondPML(m_molecule, argv[1]);

  log("rigidity")<<"Molecule has:"<<endl;
  log("rigidity") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("rigidity")<<"> "<<protein->m_initialCollisions.size()<<" initial collisions"<<endl;
  log("rigidity")<<"> "<<protein->m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("rigidity") << "> " << protein->m_spanning_tree->getNumDOFs() << " DOFs of which " << protein->m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  Configuration* conf = new Configuration(protein);
  protein->setConfiguration(conf);
  //conf->computeCycleJacobianAndNullSpace();

  log("rigidity")<<"Dimension of Jacobian: " << conf->getNullspace()->getMatrix()->size1 << " rows, ";
  log("rigidity")<< conf->getNullspace()->getMatrix()->size2<<" columns"<<endl;
  log("rigidity")<<"Dimension of kernel "<< conf->getNullspace()->getNullspaceSize()<<endl;

  int sample_id = 1;
  string out_file = out_path + "output/" + name + "_new_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".pdb";
  if(options.saveData > 0){
    ///Write PDB File for pyMol usage
    IO::writePdb(protein, out_file);
  }
  if(options.saveData > 1){
    ///save pyMol coloring script
    string pyMol=out_path + "output/" +  name + "_pyMol_" +
                 std::to_string((long long)sample_id)
                 //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                 + ".pml";
    string statFile=out_path + "output/" +  name + "_stats_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".txt";
    ///Write pyMol script
    IO::writePyMolScript(protein, out_file, pyMol);
    ///Write statistics
    IO::writeStats(protein,statFile);

    if(options.saveData > 2){
      ///save Jacobian and Nullspace to file
      string outJac=out_path + "output/" +  name + "_jac_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".txt";
      string outNull=out_path + "output/" +  name + "_nullSpace_" +
                     std::to_string((long long)sample_id)
                     //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                     + ".txt";
      ///save singular values
      string outSing=out_path + "output/" +  name + "_singVals_" +
                     std::to_string((long long)sample_id)
                     //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                     + ".txt";
      string rbFile=out_path + "output/" +  name + "_RBs_" +
                    std::to_string((long long)sample_id)
                    //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                    + ".txt";
      string covFile=out_path + "output/" +  name + "_covBonds_" +
                     std::to_string((long long)sample_id)
                     //static_cast<ostringstream*>( &(ostringstream() << sample_id) )->str()
                     + ".txt";

      ///Write Jacobian
      //gsl_matrix_outtofile(m_protein->m_conf->CycleJacobian, outJac);
      /////Write Null Space matrix
      //gsl_matrix_outtofile(conf->CycleNullSpace->m_nullspaceBasis,outNull);
      /////Write Singular values
      //gsl_vector_outtofile(NullSpaceRet::singularValues,outSing);
      if (NullspaceSVD* derived = dynamic_cast<NullspaceSVD*>(conf->getNullspace())) {
        derived->writeMatricesToFiles(outJac, outNull, outSing);
      }

      ///Write rigid bodies
      IO::writeRBs(protein, rbFile);
      ///Write covalent bonds
      IO::writeCovBonds(protein,covFile);
    }
  }

  return 0;
}


