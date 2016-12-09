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
#include <math/NullspaceSVD.h>
#include <gsl/gsl_vector_double.h>

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
#include "core/HBond.h"
#include "Logger.h"
#include "applications/options/ExploreOptions.h"
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

  //ExploreOptions options(argc,argv);
  ExploreOptions::createOptions(argc, argv);
  ExploreOptions &options = *(ExploreOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//ExploreOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      movingResidues,
      options.extraCovBonds,
      options.roots,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  protein->setCollisionFactor(options.collisionFactor);

  string out_path = options.workingDirectory;
  string name = protein->getName();

  string hBondOut = "hBonds_out.txt";
  string hBondIn = "hBonds_in.txt";
  IO::writeHbonds(protein,hBondOut );
  IO::writeHbondsIn(protein,hBondIn );

  Configuration* conf = protein->m_conf;

  double initialHbondEnergy = HbondIdentifier::computeHbondEnergy(conf);
  double initialVdwEnergy = protein->vdwEnergy(options.collisionCheck);
  //conf->computeCycleJacobianAndNullSpace();

  log("hierarchy") << "Molecule has:" << endl;
  log("hierarchy") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("hierarchy") << "> " << protein->getInitialCollisions().size() << " initial collisions" << endl;
  log("hierarchy") << "> " << protein->m_spanningTree->m_cycleAnchorEdges.size() << " hydrogen bonds" << endl;
  log("hierarchy") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " <<
  protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  NullspaceSVD ns = *(dynamic_cast<NullspaceSVD*>(conf->getNullspace()));
  //Copy to local variable
  int numCols = ns.getMatrix()->size2;
  int nullspaceCols = ns.getNullspaceSize();
  int sampleCount = 0;
  gsl_vector* singValVector = gsl_vector_copy(ns.getSVD()->S);
  gsl_matrix* baseNullspaceV = gsl_matrix_calloc(ns.getSVD()->V->size1,ns.getSVD()->V->size2);
  gsl_matrix_memcpy(baseNullspaceV, ns.getSVD()->V);
  gsl_matrix* baseJacobian = gsl_matrix_calloc(ns.getMatrix()->size1,ns.getMatrix()->size2);
  gsl_matrix_memcpy(baseJacobian, ns.getMatrix());

  log("hierarchy") << "Dimension of Jacobian: " << ns.getMatrix()->size1 << " rows, ";
  log("hierarchy") << numCols << " columns" << endl;
  log("hierarchy") << "Dimension of kernel " << nullspaceCols << endl;

  log("hierarchy") << "Initial hbond energy: " << initialHbondEnergy << endl << endl;

  if(options.saveData > 1) {
    string out_file = out_path + "output/" + name + ".pdb";
    ///save pyMol coloring script
    string pyMol = out_path + "output/" + name + "_pyMol.pml";
    string statFile = out_path + "output/" + name + "_stats.txt";
    string rbFile=out_path + "output/" +  name + "_RBs.txt";
    string singVals = out_path + "output/singVals.txt";
    ///Write pyMol script
    IO::writePyMolScript(protein, out_file, pyMol);
    ///Write statistics
    IO::writeStats(protein, statFile);
    ///Write rigid bodies
    IO::writeRBs(protein, rbFile);
    gsl_vector_outtofile(singValVector, singVals);
  }
//  cout<<"First conf "<<conf<<", S: "<<conf->getNullspace()->getSVD()->S<<endl;

  gsl_vector* projected_gradient = gsl_vector_calloc(numCols);
  gsl_vector* allDofs = gsl_vector_calloc(protein->m_spanningTree->getNumDOFs());

  //Write the complete J*V product out to file
  gsl_matrix* fullProduct = gsl_matrix_alloc(baseJacobian->size1, numCols);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, baseJacobian, baseNullspaceV, 0.0, fullProduct);
  string outProd="fullProduct_JV.txt";
  gsl_matrix_outtofile(fullProduct, outProd);
  gsl_matrix_free(fullProduct);

  string outMat="Vmatrix.txt";
  gsl_matrix_outtofile(baseNullspaceV, outMat);

  //Store output data in this file, space-separated in this order
  log("data")<<"sample inCollision inNullspace gradientNorm violation hbondDelta"<<endl;

  for( int v_i = 0; v_i < numCols; ++v_i) {
    conf->updateMolecule();
    bool inNullspace = v_i< nullspaceCols;
    if( v_i == nullspaceCols){
      log("hierarchy")<<endl<<"Now motions outside of the nullspace."<<endl<<endl;
    }

    gsl_vector_view projected_gradient_view = gsl_matrix_column(baseNullspaceV,numCols - v_i - 1);
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

    //Identify predictedViolation
//    double predictedViolation = gsl_vector_get(singValVector,numCols - v_i - 1)*options.stepSize;
    //Use the multiplication to ensure it works, even for systems with 5m < n (where the singular value does not exist)
    gsl_vector* violationVec= gsl_matrix_vector_mul(baseJacobian,projected_gradient);
    double predictedViolation = gsl_vector_length(violationVec);
    gsl_vector_free(violationVec);

    double gradNorm = gsl_vector_length(projected_gradient);
    //Now we have the correct cycle-dof projected gradient --> we need to scale it to the full-dof vector=
    // Convert back to full length DOFs vector
    for( auto const& edge: protein->m_spanningTree->Edges){
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
    qNew->m_id = v_i+1;
    std::copy(
        allDofs->data,
        allDofs->data + qNew->getNumDOFs(),
        qNew->m_dofs);

    string outFile = "output/allDofs_"+std::to_string(static_cast<long long>(v_i+1))+".txt";
    gsl_vector_outtofile(allDofs, outFile);

    bool inCollision = qNew->updatedMolecule()->inCollision();
    if (inCollision) {
      log("hierarchy") << "Configuration in m_direction " << v_i+1 << " is in collision. " << endl;
    }
//    else {//collision-free //todo: do we want to reject colliding configurations?
    qNew->updateMolecule();

    //Potentially reject new config if large violations?
    double observedViolation = protein->checkCycleClosure(qNew);

    /// New test with reclosing the cycles to obtain reduced energy violations etc.
//    double eps = 1e-12;
//    int maxIts = 10;
//    gsl_vector* currentViolation = gsl_vector_alloc(qNew->getCycleJacobian()->size1);
//    gsl_vector* qSol = gsl_vector_calloc(protein->m_spanningTree->getNumDOFs());
//    int count = 0;
//    while(observedViolation > eps && count <maxIts){
//      cout<<"Iteration "<<count++<<", violation: "<<observedViolation<<endl;
//      protein->computeCycleViolation(qNew,currentViolation);//compute residuum vector
//      NullspaceSVD* ns_svd = static_cast<NullspaceSVD*> (qNew->getNullspace() ); //updates Jacobian and nullspace
//      ns_svd->updateFromMatrix();
//      gsl_matrix* Jinv = ns_svd->getSVD()->PseudoInverse();
//      cout<<"Jinv dims: "<<Jinv->size1<<", "<<Jinv->size2<<endl;
//      //delta vetor, on small dimensions
//      gsl_vector* deltaQ = gsl_matrix_vector_mul( Jinv, currentViolation );
//
//      for( auto const& edge: protein->m_spanningTree->Edges){
//        int dof_id = edge->getDOF()->getIndex();
//        int cycle_dof_id = edge->getDOF()->getCycleIndex();
//        if ( cycle_dof_id!=-1 ) {
//          gsl_vector_set(qSol,dof_id,gsl_vector_get(allDofs,dof_id) - gsl_vector_get(deltaQ,cycle_dof_id));
//        }
//        else if ( dof_id!=-1 ) {//use zeros for free dofs
//          gsl_vector_set(qSol,dof_id,0);
//        }
//      }
//      gsl_vector_free(deltaQ);
//      for (int ind = 0; ind < qSol->size; ++ind)
//        gsl_vector_set(allDofs, ind, gsl_vector_get(qSol,ind));
//
//      delete qNew; //delete previous trial
//      Configuration *qNew = new Configuration(conf);
//      qNew->m_id = v_i+1;
//      std::copy(
//          qSol->data,
//          qSol->data + qNew->getNumDOFs(),
//          qNew->m_dofs);
//
//      qNew->updateMolecule();
//      observedViolation = protein->checkCycleClosure(qNew);
//    }
//    gsl_vector_free(currentViolation);
//    gsl_vector_free(qSol);
//    /// END NEW TEST

    qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(ExploreOptions::getOptions()->collisionCheck);
//    double hBondEnergy = HbondIdentifier::computeHbondEnergy(qNew);
//    double normDeltaHEnergy = hBondEnergy - initialHbondEnergy;
    double normDeltaHEnergy = HbondIdentifier::computeHbondNormedEnergyDifference(qNew);
    double deltaVdwEnergy = qNew->getMolecule()->vdwEnergy(options.collisionCheck) - initialVdwEnergy;

    qNew->writeQToBfactor();
    log("hierarchy") << "> New structure: " << ++sampleCount;
    log("hierarchy") << " of a total of " << ns.getMatrix()->size2 << " samples.";
    log("hierarchy") << " Delta hbond energy: " << normDeltaHEnergy;
    log("hierarchy") << ", pred. violation: "<<predictedViolation;
    log("hierarchy") << ", obs. violation: "<<observedViolation;
    log("hierarchy") << ", delta vdw: "<<deltaVdwEnergy<<endl;
    IO::writeNewSample(qNew, conf, sampleCount, options.workingDirectory, options.saveData);

    hBondOut = "output/hBonds_"+std::to_string(static_cast<long long>(v_i+1))+".txt";
    IO::writeHbondsChange(qNew,hBondOut);

    //Store output data in this file, space-separated in this order
//    log("data")<<"sample inCollision inNullspace gradientNorm predictedViolation observedViolation hbondDelta"<<endl;
    log("data")<<sampleCount<<" "<<inCollision<<" "<<inNullspace<<" "<<gradNorm<<" "<<predictedViolation<<" "<<observedViolation<<" "<<normDeltaHEnergy<<endl;
  }
  gsl_vector_free(projected_gradient);
  gsl_vector_free(allDofs);
  gsl_vector_free(singValVector);
  gsl_matrix_free(baseJacobian);
  gsl_matrix_free(baseNullspaceV);


  reportStream.close();
  dataStream.close();

  log("hierarchy")<<"Done."<<endl;

  return 0;
}
