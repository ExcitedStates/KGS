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
#include "core/Molecule.h"
#include "core/Grid.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"

#include "planners/SamplingPlanner.h"
#include "planners/RRTPlanner.h"
#include "planners/DihedralRRT.h"
#include "planners/PoissonPlanner.h"
#include "planners/BidirectionalMovingFront.h"
#include <moves/NullspaceMove.h>
#include <moves/ClashAvoidingMove.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <moves/DecreaseStepMove.h>
#include <metrics/RMSDnosuper.h>
#include <planners/PoissonPlanner2.h>
#include <applications/options/ExploreOptions.h>
#include <math/NullspaceSVD.h>
#include <planners/MCMCPlanner.h>

using namespace std;

extern double jacobianAndNullspaceTime;
extern double rigidityTime;
extern double selectNodeTime;
//extern double SINGVAL_TOL;

void randomSampling(ExploreOptions& options) {

  Selection movingResidues(options.residueNetwork);
  Molecule *protein = IO::readPdb(
      options.initialStructureFile,
      movingResidues,
      options.extraCovBonds,
      options.roots,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  protein->setCollisionFactor(options.collisionFactor);

  if (options.collapseRigid > 0) {
    log("samplingStatus") << "Before collapsing edges:" << endl;
    log("samplingStatus") << "> " << protein->m_spanningTree->m_cycleAnchorEdges.size() << " constraints" << endl;
    log("samplingStatus") << "> " << protein->m_spanningTree->getNumDOFs()  <<  " DOFs of which "  <<  protein->m_spanningTree->getNumCycleDOFs()  <<  " are cycle-DOFs"  <<  endl;
    log("samplingStatus") << "> " << protein->m_spanningTree->Vertex_map.size()  <<  " rigid bodies"  <<  endl;
    protein = protein->collapseRigidBonds(options.collapseRigid);
  }
//  if(!options.annotationFile.empty())
//    IO::readAnnotations(protein, options.annotationFile);

  for(auto const& coll: protein->getAllCollisions()){
    log("dominik") << "Ini coll: " << coll.first->getId() << " " << coll.first->getName() << " " << coll.second->getId() << coll.second->getName() << endl;
  }

  log("samplingStatus") << "Molecule has:" << endl;
  log("samplingStatus") << "> " << protein->getAtoms().size()  <<  " atoms"  <<  endl;
  log("samplingStatus") << "> " << protein->getInitialCollisions().size() << " initial collisions" << endl;
  log("samplingStatus") << "> " << protein->m_spanningTree->m_cycleAnchorEdges.size() << " constraints" << endl;
  log("samplingStatus") << "> " << protein->m_spanningTree->getNumDOFs()  <<  " DOFs of which "  <<  protein->m_spanningTree->getNumCycleDOFs()  <<  " are cycle-DOFs"  <<  endl;
  log("samplingStatus") << "> " << protein->m_spanningTree->Vertex_map.size()  <<  " rigid bodies"  <<  endl;

  if(protein->m_spanningTree->m_cycleAnchorEdges.size()==0){
    log("samplingStatus")<<"Warning: There are no constraints"<<endl;
//    log("samplingStatus")<<"Stopping because there are no hydrogen bonds"<<endl;
//    exit(-1);
  }
  //Compute energy
  protein->m_conf->m_vdwEnergy = protein->vdwEnergy(options.collisionCheck);
  log("samplingStatus")<<"> "<<protein->m_conf->m_vdwEnergy<<" kcal/mole energy"<<endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  Selection metricSelection(options.metricSelection);
  try {
    if(ExploreOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(ExploreOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(ExploreOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new ClashAvoidingMove(options.maxRotation,
                                 options.decreaseSteps,
                                 options.collisionCheck,
                                 options.projectConstraints);
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(ExploreOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  //Initialize direction
  Direction* direction;
  Selection selectionMoving(options.gradientSelection);

  try {
    if (options.gradient == 0) {
      log("samplingStatus") << "Using random direction" << endl;
      direction = new RandomDirection(selectionMoving);
    } else if (options.gradient <= 2) {
      log("samplingStatus") << "Using dihedral direction" << endl;
      direction = new DihedralDirection(selectionMoving);
    } else if (options.gradient <= 4) {
      log("samplingStatus") << "Using MSD direction" << endl;
      direction = new MSDDirection(selectionMoving, options.alignAlways);
    } else if (options.gradient <= 5) {
      log("samplingStatus") << "Using LS direction" << endl;
      direction = new LSNullspaceDirection(selectionMoving);
    } else {
      cerr << "Unknown gradient specified" << endl;
      exit(-1);
    }
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize planner
  SamplingPlanner* planner;
  if(options.planner_string=="binnedrrt"){
    log("samplingStatus")<<"Using binned RRT planner"<<endl;
    planner = new RRTPlanner(   // protein, *direction, options.explorationRadius );
        protein,
        direction,
        options.explorationRadius,
        options.samplesToGenerate,
        options.gradient,
        options.stepSize,
        options.maxRotation,
        options.scaleToRadius
    );
  }else if(options.planner_string=="dihedralrrt"){
    log("samplingStatus")<<"Using dihedral RRT planner"<<endl;
    planner = new DihedralRRT(
        protein,
        direction,
        options.samplesToGenerate,
        options.explorationRadius,
        options.maxRotation,
        options.sampleRandom
    );
  }else if(options.planner_string=="poisson"){
    log("samplingStatus")<<"Using Poisson-disk planner"<<endl;
    planner = new PoissonPlanner(
        protein,
        options.samplesToGenerate,
        options.poissonMaxRejectsBeforeClose,
        options.stepSize,
        options.gradientSelection,
        options.enableBVH
    );
  }else if(options.planner_string=="mcmc"){
    log("samplingStatus")<<"Using MCMC planner"<<endl;
    planner = new MCMCPlanner(
        protein,
        direction,
        options.samplesToGenerate,
        options.stepSize
    );
  }else{
    cerr<<"Unknown planner option specified!"<<endl;
    exit(-1);
  }
  planner->initialize(move, metric, options.workingDirectory, options.saveData);

  log() << "Total DOFs: " << protein->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs()
      << ", Max accessible DOFs: " << protein->m_spanningTree->getNumDOFs() - protein->m_spanningTree->getNumCycleDOFs() +
      protein->m_conf->getNullspace()->getNullspaceSize() << endl;fflush(stdout);

  if(options.saveData > 1){
    string out = options.workingDirectory + "output/" + protein->getName() + "_q_0.txt";
    IO::writeQ(protein, protein->m_conf, out);
  }

  log()<<"Number of rigid clusters: "<<protein->m_conf->m_numClusters;
  log()<<", biggest cluster: index "<<protein->m_conf->m_maxIndex<<" with "<<protein->m_conf->m_maxSize<<" atoms!"<<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidDihedrals() << " rigidified";
  log()<<" and " << ( protein->m_conf->getNullspace()->getNumDOFs()-
                      protein->m_conf->getNullspace()->getNumRigidDihedrals()) << " coordinated dihedrals" <<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidHBonds()<<" rigid out of "<<protein->getHBonds().size()<<" hydrogen bonds!"<<endl<<endl;

  log("samplingStatus")<<"Sampling ...\n";
  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  //Start exploring
  planner->generateSamples();

  //Print final status
  double end_time = timer.ElapsedTime();
  std::list<Configuration*>& m_samples = planner->getSamples();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianAndNullspaceTime<<" seconds\n";
  log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
  log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";


  if(options.saveData > 0){
    planner->createTrajectory();
  }

  //Clean up
  delete planner;
  delete direction;
}

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

  ofstream reportStream;
  reportStream.open("kgs_report.log");
  enableLogger("report", reportStream);

  ofstream plannerStream;
  plannerStream.open("kgs_planner.log");
  enableLogger("dominik", plannerStream);

  ofstream debugStream;
  debugStream.open("kgs_debug.log");
  enableLogger("debug", debugStream);

  ExploreOptions::createOptions(argc, argv);

  ExploreOptions &options = *(ExploreOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//ExploreOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  // Set SVD cutoff
//  SINGVAL_TOL = options.svdCutoff;
  NullspaceSVD::setSingularValueTolerance(options.svdCutoff);

  randomSampling(options);

  reportStream.close();
  plannerStream.close();

  return 0;
}
