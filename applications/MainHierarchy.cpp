//
// Created by Dominik Budday on 06.06.16.
//

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math/gsl_helpers.h>
#include <gsl/gsl_matrix_double.h>

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

int main( int argc, char* argv[] ) {

  enableLogger("hierarchy");
  enableLogger("samplingStatus");

  ofstream reportStream;
  reportStream.open("kgs_report.log");
  enableLogger("report", reportStream);

  ofstream dataStream;
  dataStream.open("hierarchy_data.txt");
  enableLogger("data", dataStream);

  if (argc < 2) {
    cerr << "Too few arguments. Please specify PDB-file in arguments" << endl;
    exit(-1);
  }

  //SamplingOptions options(argc,argv);
  SamplingOptions::createOptions(argc, argv);
  SamplingOptions &options = *(SamplingOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//SamplingOptions
    options.print();
  }


  // Set seed
  srand(options.seed);

  Molecule protein;
  protein.setCollisionFactor(options.collisionFactor);

  IO::readPdb(&protein, options.initialStructureFile, options.extraCovBonds);
//  options.setResidueNetwork(&protein);

  string out_path = options.workingDirectory;
  string name = protein.getName();

  if(options.hydrogenbondMethod=="user")
    IO::readHbonds( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="rnaview")
    IO::readHbonds_rnaview( &protein, options.hydrogenbondFile, options.annotationFile.empty() );
  else if(options.hydrogenbondMethod=="first" || options.hydrogenbondMethod=="FIRST")
    IO::readHbonds_first( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="kinari" || options.hydrogenbondMethod=="KINARI")
    IO::readHbonds_kinari( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="hbplus" || options.hydrogenbondMethod=="hbPlus")
    IO::readHbonds_hbPlus( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="vadar")
    IO::readHbonds_vadar( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="dssr")
    IO::readHbonds_dssr( &protein, options.hydrogenbondFile );
  else if(options.hydrogenbondMethod=="identify")
    HbondIdentifier::identifyHbonds(&protein);


  string hBondOut = "hBonds_out.txt";
  string hBondIn = "../hBonds_in.txt";
  IO::writeHbonds(&protein,hBondOut );
  IO::writeHbondsIn(&protein,hBondIn );

  IO::readRigidbody(&protein);
  protein.buildSpanningTree();

  Configuration *conf = new Configuration(&protein);
  protein.setConfiguration(conf);
  conf->updateMolecule();

  protein.m_initialCollisions = protein.getAllCollisions();

  double initialHbondEnergy = HbondIdentifier::computeHbondEnergy(conf);
  double initialVdwEnergy = protein.vdwEnergy(options.collisionCheck);
  //conf->computeCycleJacobianAndNullSpace();

  log("hierarchy") << "Molecule has:" << endl;
  log("hierarchy") << "> " << protein.getAtoms().size() << " atoms" << endl;
  log("hierarchy") << "> " << protein.m_initialCollisions.size() << " initial collisions" << endl;
  log("hierarchy") << "> " << protein.m_spanning_tree->CycleAnchorEdges.size() << " hydrogen bonds" << endl;
  log("hierarchy") << "> " << protein.m_spanning_tree->getNumDOFs() << " DOFs of which " <<
  protein.m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  int numCols = conf->getNullspace()->Matrix()->size2;
  int nullspaceCols = conf->getNullspace()->NullspaceSize();
  int sampleCount = 0;

  log("hierarchy") << "Dimension of Jacobian: " << conf->getNullspace()->Matrix()->size1 << " rows, ";
  log("hierarchy") << numCols << " columns" << endl;
  log("hierarchy") << "Dimension of kernel " << nullspaceCols << endl;

  log("hierarchy") << "Initial hbond energy: " << initialHbondEnergy << endl << endl;

  if(options.saveData > 2) {
    string out_file = out_path + "output/" + name + ".pdb";
    ///save pyMol coloring script
    string pyMol = out_path + "output/" + name + "_pyMol.pml";
    string statFile = out_path + "output/" + name + "_stats.txt";
    string rbFile=out_path + "output/" +  name + "_RBs.txt";
    string singVals = out_path + "output/singVals.txt";
    ///Write pyMol script
    IO::writePyMolScript(&protein, out_file, pyMol);
    ///Write statistics
    IO::writeStats(&protein, statFile);
    ///Write rigid bodies
    IO::writeRBs(&protein, rbFile);
    gsl_vector_outtofile(conf->getNullspace()->getSVD()->S,singVals);
  }
//  cout<<"First conf "<<conf<<", S: "<<conf->getNullspace()->getSVD()->S<<endl;

  gsl_vector* projected_gradient = gsl_vector_calloc(numCols);
  gsl_vector* allDofs = gsl_vector_calloc(protein.m_spanning_tree->getNumDOFs());

  //Write the complete J*V product out to file
  gsl_matrix* fullProduct = gsl_matrix_alloc(conf->getNullspace()->Matrix()->size1, numCols);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, conf->getNullspace()->Matrix(), conf->getNullspace()->getSVD()->V, 0.0, fullProduct);
  string outProd="fullProduct_JV.txt";
  gsl_matrix_outtofile(fullProduct, outProd);

  //Store output data in this file, space-separated in this order
  log("data")<<"sample inCollision inNullspace gradientNorm violation hbondDelta"<<endl;

  for( int i = 0; i < numCols; ++i) {
    conf->updateMolecule();
    bool inNullspace = i< nullspaceCols;
    if( i == nullspaceCols){
      log("hierarchy")<<endl<<"Now motions outside of the nullspace."<<endl<<endl;
    }

    gsl_vector_view projected_gradient_view = gsl_matrix_column(conf->getNullspace()->getSVD()->V,numCols - i - 1);
    gsl_vector_memcpy(projected_gradient, &projected_gradient_view.vector);

    //Scale to desired step size
    gsl_vector_scale_to_length(projected_gradient, options.stepSize);

    //    // Control the max amount of rotation
//    double max_rotation = 0;
//    for (int j=0; j<projected_gradient->size; ++j) {
//      double abs_value = Math::Abs(gsl_vector_get(projected_gradient,j));
//      if ( abs_value > max_rotation ) {
//        max_rotation = abs_value;
//      }
//    }
//    if ( max_rotation > options.maxRotation ){//MAX_ROTATION ) {
//      gsl_vector_scale(projected_gradient, options.maxRotation/max_rotation); // MAX_ROTATION/max_rotation);
//      log("hierarchy")<<"Scaled projected gradient down to: "<<gsl_vector_length(projected_gradient)<<endl;
//    }

    //Identify correct range
//    for (int j=0; j<projected_gradient->size; ++j) {
//      gsl_vector_set(projected_gradient,j,formatRangeRadian(gsl_vector_get(projected_gradient,j)));
//    }

    //Identify predictedViolation from product with constraint Jacobian
    gsl_vector* violationVec = gsl_matrix_vector_mul(conf->getNullspace()->getSVD()->matrix, projected_gradient);
    double predictedViolation = gsl_vector_length(violationVec);
//    cout<<"Pred norm: "<<predictedViolation<<", sing vec norm: "<<gsl_vector_get(conf->getNullspace()->getSVD()->S,numCols - i - 1)<<endl;
//    cout<<"Second conf "<<conf<<", S: "<<conf->getNullspace()->getSVD()->S<<endl;
    gsl_vector_free(violationVec);

    double gradNorm = gsl_vector_length(projected_gradient);

    //Now we have the correct cycle-dof projected gradient --> we need to scale it to the full-dof vector=
    // Convert back to full length DOFs vector
    for( auto const& edge: protein.m_spanning_tree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(allDofs,dof_id,gsl_vector_get(projected_gradient,cycle_dof_id));
      }
      else if ( dof_id!=-1 ) {//use zeros for free dofs
        gsl_vector_set(allDofs,dof_id,0);
      }
    }

    Configuration *qNew = new Configuration(conf);
    qNew->m_id = i+1;
    std::copy(
        allDofs->data,
        allDofs->data + qNew->getNumDOFs(),
        qNew->m_dofs);

    string outFile = "output/allDofs_"+std::to_string(static_cast<long long>(i+1))+".txt";
    gsl_vector_outtofile(allDofs, outFile);

    bool inCollision = qNew->updatedMolecule()->inCollision();
    if (inCollision) {
      log("hierarchy") << "Configuration in direction " << i+1 << " is in collision. " << endl;
    }
//    else {//collision-free //todo: do we want to reject colliding configurations?
    qNew->updateMolecule();

    //Potentially reject new config if large violations?
    double observedViolation = protein.checkCycleClosure(qNew);
    qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(SamplingOptions::getOptions()->collisionCheck);
//    double hBondEnergy = HbondIdentifier::computeHbondEnergy(qNew);
//    double normDeltaHEnergy = hBondEnergy - initialHbondEnergy;
    double normDeltaHEnergy = HbondIdentifier::computeHbondNormedEnergyDifference(qNew);
    double deltaVdwEnergy = qNew->getMolecule()->vdwEnergy(options.collisionCheck) - initialVdwEnergy;

    log("hierarchy") << "> New structure: " << ++sampleCount << " of a total of " <<
    conf->getNullspace()->Matrix()->size2 << " samples. Delta hbond energy: " << normDeltaHEnergy<<", pred. violation: "<<predictedViolation<<", obs. violation: "<<observedViolation<<", delta vdw: "<<deltaVdwEnergy<<endl;
    SamplingPlanner::writeNewSample(qNew, conf, sampleCount);

    hBondOut = "output/hBonds_"+std::to_string(static_cast<long long>(i+1))+".txt";
    IO::writeHbondsChange(qNew,hBondOut);

    //Store output data in this file, space-separated in this order
//    log("data")<<"sample inCollision inNullspace gradientNorm predictedViolation observedViolation hbondDelta"<<endl;
    log("data")<<sampleCount<<" "<<inCollision<<" "<<inNullspace<<" "<<gradNorm<<" "<<predictedViolation<<" "<<observedViolation<<" "<<normDeltaHEnergy<<endl;
  }
  gsl_vector_free(projected_gradient);
  gsl_vector_free(allDofs);

  reportStream.close();
  dataStream.close();

  log("hierarchy")<<"Done."<<endl;

  return 0;
}
