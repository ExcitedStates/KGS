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
#include <gsl/gsl_vector.h>

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
#include "applications/options/HierarchyOptions.h"
#include "moves/CompositeMove.h"
#include "directions/MSDDirection.h"
#include "directions/LSNullspaceDirection.h"

extern double jacobianAndNullspaceTime;
extern double rigidityTime;

using namespace std;

void imodeComparisonFiles(Molecule* protein, HierarchyOptions& options, const NullspaceSVD &ns, gsl_matrix* baseNullspaceV, gsl_vector* singValVector, int numResis){
  //Save the output for comparison to normal modes via iMode
  //This requires a contiguous residue sequence without gaps
  Residue* const firstResi = (protein->getAtoms())[0]->getResidue();
  int numCols = ns.getMatrix()->size2;
  int dofCounter = 0;
  int numProlines = 0;
  gsl_vector* projected_gradient = gsl_vector_alloc(numCols);
  std::map< pair<Residue*,int>, int> iModDofIndices;
  std::vector<Residue *>::iterator resIt;
  for(auto chain : protein->m_chains) {
    std::vector<Residue *> &residues = chain->getResidues();
    for(auto const resi : residues){
      if (resi == firstResi){//Special case where first dof is fixed in iMod
        iModDofIndices.insert(make_pair( make_pair(resi,1), dofCounter++) );
        continue;
      }
      iModDofIndices.insert(make_pair( make_pair(resi,0), dofCounter++) );
      if (resi->getName()!="PRO"){
        iModDofIndices.insert(make_pair( make_pair(resi,1), dofCounter++) );
      }
      else
        numProlines++;
    }
  }
  int numImodDofs = numResis * 2 - 1 - numProlines;
  log("hierarchy")<<"KGS predicts "<<numImodDofs<<" imod dofs."<<endl;
  gsl_vector* iModDofs = gsl_vector_alloc(numImodDofs);

  //Open the output file
  string outIMode = "singVecsKGS_ic.evec";
  ofstream output( outIMode.c_str() );
  if(!output.is_open()) {
    cerr<<"Cannot write to "<<outIMode<<". You might need to create output directory first"<<endl;
    exit(-1);
  }
  //Write header
  output<<"Eigenvector file: COVAR"<<endl;
  output<<numCols<<" "<<numImodDofs<<" Contains "<<numCols<<" singular vectors"<<endl;
  //Loop through the dofs in V and in the molecule
  for( int v_i = 0; v_i < numCols; ++v_i) {
    gsl_vector_view projected_gradient_view = gsl_matrix_column(baseNullspaceV,numCols - v_i - 1);
    gsl_vector_memcpy(projected_gradient, &projected_gradient_view.vector);
    //Scale to desired step size
    gsl_vector_scale_to_length(projected_gradient, options.stepSize);
    //Refill iModDofs vector
    gsl_vector_set_zero(iModDofs);
    for( auto const edge : protein->m_spanningTree->m_edges){
      Bond *bond = edge->getBond();
      // Check for global dof bonds
      if (bond == nullptr)
        continue;

      int idx=-1; //iMod index
      int cycleDofID =  edge->getDOF()->getCycleIndex(); //KGS cycle index
      if (cycleDofID != -1) {//if cycle DOF in KGS
        Atom *atom1 = bond->m_atom1;
        if (bond->m_atom1->getName() == "N" and bond->m_atom2->getName() == "CA") {
          idx = iModDofIndices[make_pair(atom1->getResidue(),0)];
//          log("hierarchy") << "First dof in resi " << atom1->getResidue()->getId() << ": " << idx << endl;
        } else {
          if (bond->m_atom1->getName() == "CA" and bond->m_atom2->getName() == "C") {
            idx = iModDofIndices[make_pair(atom1->getResidue(),1)];
//            log("hierarchy") << "Second dof in resi " << atom1->getResidue()->getId() << ": " << idx << endl;
          } else
            idx = -1;
        }
        if (idx >= 0 && idx <= numImodDofs) {
//          log("hierarchy") << "IDX " << idx << " length " << iModDofs->size << " cycle IDX " << cycleDofID << " length "
//               << projected_gradient->size << endl;
          gsl_vector_set(iModDofs, idx, gsl_vector_get(projected_gradient, cycleDofID));
        }
      }
    }
    //Scale to desired step size
    gsl_vector_scale_to_length(iModDofs, options.stepSize);
    //Write iMod vector information to file
    output<<"****"<<endl;
    if (numCols - v_i - 1 < singValVector->size) {
      output << v_i + 1 << " " << gsl_vector_get(singValVector, numCols - v_i - 1) << endl;
    }
    else{
      output << v_i + 1 << " " << 0.0 << endl;
    }
    for(int i=0; i<numImodDofs; i++){
      output<<gsl_vector_get(iModDofs,i);
      if (i<numImodDofs-1) {
        output << " ";
      }
    }
    output<<endl;
  }
  output.close();
  gsl_vector_free(projected_gradient);
  gsl_vector_free(iModDofs);
}


int main( int argc, char* argv[] ) {

  enableLogger("hierarchy");
  enableLogger("samplingStatus");

//  ofstream reportStream;
//  reportStream.open("kgs_report.log");
//  enableLogger("report", reportStream);

  ofstream dataStream;
  dataStream.open("hierarchy_data.txt");
  enableLogger("data", dataStream);

//  ofstream debugStream;
//  debugStream.open("kgs_debug.log");
//  enableLogger("debug", debugStream);

  if (argc < 2) {
    cerr << "Too few arguments. Please specify PDB-file in arguments" << endl;
    exit(-1);
  }

  //HierarchyOptions options(argc,argv);
  HierarchyOptions::createOptions(argc, argv);
  HierarchyOptions &options = *(HierarchyOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//HierarchyOptions
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

  int numResis = 0;
  for(auto chain : protein->m_chains) {
    std::vector<Residue *> &residues = chain->getResidues();
    numResis += residues.size();
  }
  log("hierarchy") << "Molecule has:" << endl;
  log("hierarchy") << "> " << protein->m_chains.size() << " chains" << endl;
  log("hierarchy") << "> " << numResis << " residues" << endl;
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

//  gsl_matrix_outtofile(baseJacobian, "gasJacobian.txt");
//  gsl_matrix_outtofile(baseNullspaceV,"gasV.txt");

  log("hierarchy") << "Dimension of Jacobian: " << ns.getMatrix()->size1 << " rows, ";
  log("hierarchy") << numCols << " columns" << endl;
  log("hierarchy") << "Dimension of kernel " << nullspaceCols << endl;
  log("hierarchy") << "Initial hbond energy: " << initialHbondEnergy << endl;
  log("hierarchy") << "Step size: " << options.stepSize << endl << endl;

  string out_file = out_path + "output/" + name + ".pdb";

  protein->writeRigidbodyIDToBFactor();
  IO::writePdb(protein,out_file);
  //Necessary for postprocessing, potentially do rigid cluster decomp as well or plot violations
  string statFile = out_path + "output/" + name + "_stats.txt";
  string singVals = out_path + "output/singVals.txt";
  IO::writeStats(protein, statFile);
  gsl_vector_outtofile(singValVector, singVals);

//  cout<<"First conf "<<conf<<", S: "<<conf->getNullspace()->getSVD()->S<<endl;

  gsl_vector* projected_gradient = gsl_vector_calloc(numCols);
  gsl_vector* allDofs = gsl_vector_calloc(protein->m_spanningTree->getNumDOFs());

  //Write the complete J*V product out to file
  gsl_matrix* fullProduct = gsl_matrix_alloc(baseJacobian->size1, numCols);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, baseJacobian, baseNullspaceV, 0.0, fullProduct);
  if(options.saveData > -1) {
    string outMat = "Vmatrix.txt";
    gsl_matrix_outtofile(baseNullspaceV, outMat);
    string outProd = "fullProduct_JV.txt";
    gsl_matrix_outtofile(fullProduct, outProd);
  }

  if(options.saveData > 1) {
    imodeComparisonFiles(protein, options, ns, baseNullspaceV, singValVector, numResis);
  }

  if(options.saveData > 2) {
    string outJac = out_path + "output/" + name + "_jac.txt";
    gsl_matrix_outtofile(baseJacobian, outJac);
    ///save pyMol coloring script
    string pyMol = out_path + "output/" + name + "_pyMol.pml";
    string rbFile=out_path + "output/" +  name + "_RBs.txt";
    ///Write pyMol script
    IO::writePyMolScript(protein, out_file, pyMol);
    ///Write statistics
    ///Write rigid bodies
    IO::writeRBs(protein, rbFile);
  }

  gsl_matrix_free(fullProduct);
  //Store output data in this file, space-separated in this order
  log("data")<<"sample inCollision inNullspace gradientNorm predictedViolation observedViolation hbondDelta"<<endl;

  int maxSamples = min(options.samples,numCols);

  for( int v_i = 0; v_i < maxSamples; ++v_i) {
    conf->updateMolecule();
    bool inNullspace = v_i< nullspaceCols;
    if( v_i == nullspaceCols){
      log("hierarchy")<<endl<<"Now motions outside of the nullspace."<<endl<<endl;
    }

    gsl_vector_view projected_gradient_view = gsl_matrix_column(baseNullspaceV,numCols - v_i - 1);
    gsl_vector_memcpy(projected_gradient, &projected_gradient_view.vector);

    //Scale to desired step size
    gsl_vector_scale_to_length(projected_gradient, options.stepSize);

    //Identify predictedViolation
//    double predictedViolation = gsl_vector_get(singValVector,numCols - v_i - 1)*options.stepSize;
    //Use the multiplication to ensure it works, even for systems with 5m < n (where the singular value does not exist)
    gsl_vector* violationVec= gsl_matrix_vector_mul(baseJacobian,projected_gradient);
    double predictedViolation = gsl_vector_length(violationVec);
    gsl_vector_free(violationVec);

    double gradNorm = gsl_vector_length(projected_gradient);
    //Now we have the correct cycle-dof projected gradient --> we need to scale it to the full-dof vector=
    // Convert back to full length DOFs vector
    for( auto const& edge: protein->m_spanningTree->m_edges){
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

    if(options.saveData > 1) {
      string outFile = "output/allDofs_" + std::to_string(static_cast<long long>(v_i + 1)) + ".txt";
      gsl_vector_outtofile(allDofs, outFile);
    }

    bool inCollision = qNew->updatedMolecule()->inCollision();
    if (inCollision) {
      log("hierarchy") << "Configuration in m_direction " << v_i+1 << " is in collision. " << endl;
    }
//    else {//collision-free //todo: do we want to reject colliding configurations?
    qNew->updateMolecule();

    //Potentially reject new config if large violations?
    double observedViolation = protein->checkCycleClosure(qNew);

    qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(HierarchyOptions::getOptions()->collisionCheck);
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

    if(options.saveData > 1) {
      hBondOut = "output/hBonds_" + std::to_string(static_cast<long long>(v_i + 1)) + ".txt";
      IO::writeHbondsChange(qNew, hBondOut);
    }

    //Store output data in this file, space-separated in this order
//    log("data")<<"sample inCollision inNullspace gradientNorm predictedViolation observedViolation hbondDelta"<<endl;
    log("data")<<sampleCount<<" "<<inCollision<<" "<<inNullspace<<" "<<gradNorm<<" "<<predictedViolation<<" "<<observedViolation<<" "<<normDeltaHEnergy<<endl;
  }

  ///If target structure is available, also test direction towards target
  if(options.targetStructureFile != ""){
    Molecule* target = IO::readPdb(
        options.targetStructureFile,
        movingResidues,
        options.extraCovBonds,
        options.roots,
        options.hydrogenbondMethod,
        options.hydrogenbondFile,
        protein
    );
    target->setCollisionFactor(options.collisionFactor);

    string targetName = target->getName();
    IO::writeHbonds(target,"hBonds_target.txt" );

    Configuration* targetConf = target->m_conf;
    Selection gradientSelection("heavy"); //maybe backbone is useful for amide-only HDXMS data comparison?
    Direction* targetGrad;
    //exchange this for other directions if desired
    targetGrad = new MSDDirection(gradientSelection, true);
    //Compute full-length gradient
    targetGrad->gradient(conf,targetConf,allDofs);

    conf->convertAllDofsToCycleDofs(projected_gradient, allDofs);

    //Scale to desired step size
    gsl_vector_scale_to_length(projected_gradient, options.stepSize);
    gsl_vector_outtofile(projected_gradient,"desiredGradient.txt");

    //Identify predictedViolation
//    double predictedViolation = gsl_vector_get(singValVector,numCols - v_i - 1)*options.stepSize;
    //Use the multiplication to ensure it works, even for systems with 5m < n (where the singular value does not exist)
    gsl_vector* violationVec= gsl_matrix_vector_mul(baseJacobian,projected_gradient);
    double predictedViolation = gsl_vector_length(violationVec);
    gsl_vector_outtofile(violationVec, "gradientVectorViolation.txt");
    gsl_vector_free(violationVec);

    //Overlap with singular vectors
    gsl_vector* overlapVec= gsl_matrix_vector_mul(baseNullspaceV,projected_gradient);
    double overlapNorm = gsl_vector_length(overlapVec);
    gsl_vector_outtofile(overlapVec, "gradientOverlap.txt");
    gsl_vector_free(overlapVec);

    double gradNorm = gsl_vector_length(projected_gradient);
    //Now we have the correct cycle-dof projected gradient --> we need to scale it to the full-dof vector=
    // Convert back to full length DOFs vector
    for( auto const& edge: protein->m_spanningTree->m_edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(allDofs,dof_id,gsl_vector_get(projected_gradient,cycle_dof_id));
      }
      else if ( dof_id!=-1 ) {//use zeros for free dofs, they don't affect h-bonds
        gsl_vector_set(allDofs,dof_id,0);
      }
    }

    Configuration *qNew = new Configuration(conf);
    qNew->m_id = maxSamples+1; //greater than previous sample
    std::copy(
        allDofs->data,
        allDofs->data + qNew->getNumDOFs(),
        qNew->m_dofs);

    string outFile = "output/targetGradDofs_"+std::to_string(static_cast<long long>(qNew->m_id))+".txt";
    gsl_vector_outtofile(allDofs, outFile);

    bool inCollision = qNew->updatedMolecule()->inCollision();
    if (inCollision) {
      log("hierarchy") << "Configuration in m_direction " << qNew->m_id << " is in collision. " << endl;
    }
//    else {//collision-free //todo: do we want to reject colliding configurations?
    qNew->updateMolecule();

    //Potentially reject new config if large violations?
    double observedViolation = protein->checkCycleClosure(qNew);
    double normDeltaHEnergy = HbondIdentifier::computeHbondNormedEnergyDifference(qNew);

    qNew->writeQToBfactor();
    log("hierarchy") << "> Target-gradient structure: " << ++sampleCount;
    log("hierarchy") << ", pred. violation: "<<predictedViolation;
    log("hierarchy") << ", obs. violation: "<<observedViolation;
    IO::writeNewSample(qNew, conf, sampleCount, options.workingDirectory, options.saveData);

    hBondOut = "output/hBonds_target_"+std::to_string(static_cast<long long>(qNew->m_id))+".txt";
    IO::writeHbondsChange(qNew,hBondOut);

    //Store output data in this file, space-separated in this order
//    log("data")<<"sample inCollision inNullspace gradientNorm predictedViolation observedViolation hbondDelta"<<endl;
    log("data")<<sampleCount<<" "<<inCollision<<" 0 "<<gradNorm<<" "<<predictedViolation<<" "<<observedViolation<<" "<<normDeltaHEnergy<<endl;
  }
  gsl_vector_free(projected_gradient);
  gsl_vector_free(allDofs);
  gsl_vector_free(singValVector);
  gsl_matrix_free(baseJacobian);
  gsl_matrix_free(baseNullspaceV);

//  reportStream.close();
  dataStream.close();

  log("hierarchy")<<"Done."<<endl;

  return 0;
}
