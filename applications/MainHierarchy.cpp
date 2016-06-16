//
// Created by Dominik Budday on 06.06.16.
//

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#include <gsl/gsl_matrix.h>
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

  string out_path = options.workingDirectory;
  //string pdb_file = path + protein_name + ".pdb";

  Molecule *protein = new Molecule();
  protein->setCollisionFactor(options.collisionFactor);

  IO::readPdb(protein, options.initialStructureFile, options.extraCovBonds);
  string name = protein->getName();

  if (options.hydrogenbondMethod == "user")
    IO::readHbonds(protein, options.hydrogenbondFile);
  else
    HbondIdentifier::identifyHbonds(protein);

  IO::readRigidbody(protein);
  protein->buildSpanningTree();
//	if(options.hydrogenbondMethod!="user")
//		writeHBondPML(m_molecule, argv[1]);

  Configuration *conf = new Configuration(protein);
  protein->setConfiguration(conf);

  protein->m_initialCollisions = protein->getAllCollisions();

  double initialHbondEnergy = HbondIdentifier::computeHbondEnergy(conf);
  //conf->computeCycleJacobianAndNullSpace();

  log("hierarchy") << "Molecule has:" << endl;
  log("hierarchy") << "> " << protein->atoms.size() << " atoms" << endl;
  log("hierarchy") << "> " << protein->m_initialCollisions.size() << " initial collisions" << endl;
  log("hierarchy") << "> " << protein->m_spanning_tree->CycleAnchorEdges.size() << " hydrogen bonds" << endl;
  log("hierarchy") << "> " << protein->m_spanning_tree->getNumDOFs() << " DOFs of which " <<
  protein->m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  int numCols = conf->getNullspace()->Matrix()->size2;
  int nullspaceCols = conf->getNullspace()->NullspaceSize();
  int sampleCount = 0;

  log("hierarchy") << "Dimension of Jacobian: " << conf->getNullspace()->Matrix()->size1 << " rows, ";
  log("hierarchy") << numCols << " columns" << endl;
  log("hierarchy") << "Dimension of kernel " << nullspaceCols << endl;

  log("hierarchy") << "Initial hbond energy: " << initialHbondEnergy << endl << endl;

  gsl_vector* projected_gradient = gsl_vector_calloc(numCols);

  for( int i = 0; i < numCols; ++i) {

    if( i == nullspaceCols){
      log("hierarchy")<<endl<<"Now motions outside of the nullspace."<<endl<<endl;
    }

    gsl_vector_view projected_gradient_view = gsl_matrix_column(conf->getNullspace()->getSVD()->V,numCols - i - 1);
    gsl_vector_memcpy(projected_gradient, &projected_gradient_view.vector);

    gsl_vector* violation = gsl_matrix_vector_mul(conf->getNullspace()->getSVD()->matrix, projected_gradient);
    log("hierarchy")<<"Violation from column "<<i+1<<": "<<gsl_vector_length(violation)<<endl;
    gsl_vector_free(violation);

    gsl_vector_scale_to_length(projected_gradient, options.stepSize);
    // Control the max amount of rotation
    for (int j=0; j<projected_gradient->size; ++j) {
      gsl_vector_set(projected_gradient,j,formatRangeRadian(gsl_vector_get(projected_gradient,j)));
//		log("dominik")<<"Projected gradient at "<<i<<": "<<gsl_vector_get(projected_gradient, i)<<endl;
    }
    double max_rotation = 0;
    for (int j=0; j<projected_gradient->size; ++j) {
      double abs_value = Math::Abs(gsl_vector_get(projected_gradient,j));
      if ( abs_value > max_rotation ) {
        max_rotation = abs_value;
      }
    }
    if ( max_rotation > options.maxRotation ){//MAX_ROTATION ) {
      gsl_vector_scale(projected_gradient, options.maxRotation/max_rotation); // MAX_ROTATION/max_rotation);
      log("hierarchy")<<"Scaled projected gradient down to: "<<gsl_vector_length(projected_gradient)<<endl;
    }

    log("hierarchy")<<" Gradient length: "<<gsl_vector_length(projected_gradient)<<endl;

    Configuration *qNew = new Configuration(conf);
    qNew->m_id = i+1;
    std::copy(
        projected_gradient->data,
        projected_gradient->data + qNew->getNumDOFs(),
        qNew->m_dofs);

    if (qNew->updatedMolecule()->inCollision()) {
      log("hierarchy") << "Configuration in direction " << i+1 << " is in collision. " << endl;
    }
//    else {//collision-free //todo: do we want to reject colliding configurations?
    protein = qNew->updatedMolecule();

    //Potentially reject new config if large violations?
    protein->checkCycleClosure(qNew);
    qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(SamplingOptions::getOptions()->collisionCheck);
    double hBondEnergy = HbondIdentifier::computeHbondEnergy(qNew);

    log("hierarchy") << "> New structure: " << ++sampleCount << " of a total of " <<
    conf->getNullspace()->Matrix()->size2 << " samples. Delta hbond energy: " << initialHbondEnergy-hBondEnergy<< endl;
    SamplingPlanner::writeNewSample(qNew, conf, sampleCount);
  }
  gsl_vector_free(projected_gradient);

  log("hierarchy")<<"Done."<<endl;

  return 0;
}