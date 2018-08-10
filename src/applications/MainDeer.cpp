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
// Created by Dominik Budday on 02.08.18.
//

#include <stdexcept>
#include <string>
#include <iostream>
#include <list>
#include <math/NullspaceSVD.h>
#include "math/gsl_helpers.h"
#include <moves/ClashAvoidingMove.h>
#include "moves/NullspaceMove.h"
#include "moves/DecreaseStepMove.h"
#include "moves/RawMove.h"
#include "moves/LSNclashAvoidingMove.h"
#include <planners/SamplingPlanner.h>
#include <planners/PoissonPlanner.h>
#include <planners/BidirectionalMovingFront.h>
#include <planners/RRTPlanner.h>
#include <planners/DihedralRRT.h>
#include <planners/DEERPlanner.h>
#include <directions/AtomPairDistanceDirection.h>
#include <directions/LSNrelativeDirection.h>
#include <directions/RelativeMSDDirection.h>
#include "core/Grid.h"
#include "core/DBond.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "directions/RandomDirection.h"
#include "directions/DihedralDirection.h"
#include "directions/MSDDirection.h"
#include "directions/LSNullspaceDirection.h"
#include "directions/BlendedDirection.h"
#include "metrics/RMSDnosuper.h"
#include "applications/options/DeerOptions.h"

using namespace std;

extern double jacobianAndNullspaceTime;
extern double rigidityTime;
extern double selectNodeTime;

void scale_gradient(gsl_vector* gradient, double maxRotation);
double dist_to_objective(std::vector< std::tuple<Atom*, Atom*, double> > &goal_distances);
void printStrain(const Molecule& mol, const DeerOptions& options);

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

  ofstream reportStream;
  reportStream.open("kgs_report.log");
  enableLogger("report", reportStream);

  ofstream plannerStream;
  plannerStream.open("kgs_planner.log");
  enableLogger("planner", plannerStream);

  ofstream debugStream;
  debugStream.open("kgs_debug.log");
  enableLogger("debug", debugStream);

  DeerOptions options(argc, argv);

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//DeerOptions
    options.print();
  }
  // Set seed
  srand(options.seed);
  // Set SVD cutoff
  NullspaceSVD::setSingularValueTolerance(options.svdCutoff);

  Selection resNetwork(options.residueNetwork);
  Molecule* protein = IO::readPdb( options.initialStructureFile, options.extraCovBonds );

  ///Now we have all constraints and desired roots present to build the tree
  protein->initializeTree(resNetwork,options.collisionFactor,options.roots);

  // Check for collision
  // This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling
  for(auto const& coll: protein->getInitialCollisions()){
    log("planner")<<"Ini coll: "<<coll.first->getId()<<" "<<coll.first->getName()<<" "<<coll.second->getId()<<coll.second->getName()<<endl;
  }

  Molecule* origProtein = protein;
  if(options.collapseRigid>0) {
    log("samplingStatus")<<"Before collapsing"<<endl;
    log("samplingStatus")<<"Molecule has:"<<endl;
    log("samplingStatus")<<"> " << protein->getAtoms().size() << " atoms" << endl;
    log("samplingStatus")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" total bond constraints"<<endl;
    log("samplingStatus")<<"> "<<protein->getHBonds().size()<<" hydrogen bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->getHydrophobicBonds().size()<<" hydrophobic bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->getDBonds().size()<<" distance bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;
    gsl_matrix_outtofile(protein->m_conf->getCycleJacobian(), "nonCollapsedCycleJacobian.txt");
    gsl_matrix_outtofile(protein->m_conf->getNullspace()->getBasis(),"nonCollapsedNullspace.txt");

    protein = protein->collapseRigidBonds(options.collapseRigid);
    protein->writeRigidbodyIDToBFactor();
    string fname = "output/collapsed.pdb";
    string fnamePml = "output/collapsed.pml";
    IO::writePdb(protein, fname);
    IO::writePyMolScript(protein,fname,fnamePml);
  }

  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" total bond constraints"<<endl;

  log("samplingStatus")<<"> "<<protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  if(options.predictStrain){
    printStrain(*protein, options);
  }

  //Initialize metric
  metrics::Metric* metric = nullptr;
  Selection metricSelection(options.metricSelection);
  try {
    if(options.metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(options.metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(options.metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

//  Initialize direction
  Selection blendedSelection("all");
  Selection gradientSelection(options.gradientSelection);


//Initialize directed DEER direction
  std::vector< std::tuple<Atom*, Atom*, double> > goal_distances =
      IO::readRelativeDistances(options.relativeDistances, protein);
  LSNrelativeDirection* direction = new LSNrelativeDirection(gradientSelection, goal_distances);

  //Initialize move
  Move* move;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding LSN move"<<endl;
    move = new LSNclashAvoidingMove(direction,
                                    options.maxRotation,
                                    options.decreaseSteps,
                                    options.collisionCheck,
                                    options.projectConstraints);
  }else{
    log("samplingStatus")<<"Using raw move with LS direction"<<endl;
    move = new RawMove(options.maxRotation);
     if(options.decreaseSteps>0){
       log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
       move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
     }
  }

  ///Initialize the randomized portion of the planner
  Direction* randDirection;
  bool blendedDir = false;

//  DihedralDirection* randDirection = new DihedralDirection(gradientSelection);
  
  if(options.gradient == 0)
    randDirection = new RandomDirection(gradientSelection);
  else if(options.gradient == 1)//this one is desired and the default
    randDirection = new DihedralDirection(gradientSelection);
  else if(options.gradient == 2){
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(new DihedralDirection(gradientSelection),0);
    bdir->addDirection(new RandomDirection(blendedSelection,options.maxRotation), 1);
    randDirection = bdir;
    blendedDir = true;
  }
  else if(options.gradient == 3)
    randDirection = new MSDDirection(gradientSelection);
  else if(options.gradient == 4){
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(new MSDDirection(gradientSelection),0);
    bdir->addDirection(new RandomDirection(blendedSelection,options.maxRotation), 1);
    randDirection = bdir;
    blendedDir = true;
  }
  else if(options.gradient <= 5) {
    randDirection = new LSNullspaceDirection(gradientSelection);
  }
  else if(options.gradient <= 6) {
    std::vector< std::tuple<Atom*, Atom*, double> > goal_distances =
        IO::readRelativeDistances(options.relativeDistances, protein);
    randDirection = new LSNrelativeDirection(gradientSelection, goal_distances);
  }
  else if(options.gradient <= 7) {
    std::vector< std::tuple<Atom*, Atom*, double> > goal_distances =
        IO::readRelativeDistances(options.relativeDistances, protein);
    Direction* d1 = new LSNrelativeDirection(gradientSelection, goal_distances);
    Direction* d2 = new RandomDirection(blendedSelection, options.maxRotation);
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(d1, 0);
    bdir->addDirection(d2, 1);
    randDirection = bdir;
    blendedDir = true;
  }
  else if(options.gradient <= 8) {
    std::vector< std::tuple<Atom*, Atom*, double> > goal_distances =
        IO::readRelativeDistances(options.relativeDistances, protein);
    Direction* d1 = new RelativeMSDDirection(goal_distances);
    Direction* d2 = new RandomDirection(blendedSelection, options.maxRotation);
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(d1, 0);
    bdir->addDirection(d2, 1);
    randDirection = bdir;
    blendedDir = true;
  }

  ///Randomized move
  Move* randMove;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding nullspace move for randomized portion"<<endl;
    randMove = new ClashAvoidingMove(options.maxRotation,
                                     options.decreaseSteps,
                                     options.collisionCheck,
                                     options.projectConstraints);
  }else{
    log("samplingStatus")<<"Using regular nullspace move for randomized portion"<<endl;
    randMove = new NullspaceMove(options.maxRotation);
    //Todo: reenable step-size if desired, currently not functional
     if(options.decreaseSteps>0){
       log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
       randMove = new DecreaseStepMove(randMove, (unsigned int)options.decreaseSteps, options.decreaseFactor);
     }
  }

  log() << "Total DOFs: " << protein->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs() << endl;fflush(stdout);
  log()<<"Number of rigid clusters: "<<protein->m_conf->m_numClusters;
  log()<<", biggest cluster: index "<<protein->m_conf->m_maxIndex<<" with "<<protein->m_conf->m_maxSize<<" atoms!"<<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidDihedrals() << " rigidified";
  log()<<" and " << ( protein->m_conf->getNullspace()->getNumDOFs()-
                      protein->m_conf->getNullspace()->getNumRigidDihedrals()) << " coordinated dihedrals" <<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidHBonds()<<" rigid out of "<<protein->getHBonds().size()<<" hydrogen bonds!"<<endl<<endl;

  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

/*  //%%%%%%%%%%%%%%%%%%%%%%DOWNHILL TO TARGET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  log("samplingStatus")<<"Initial distance to objective: "<<dist_to_objective(goal_distances)<<endl;
//  log("samplingStatus")<<"Performing downhill motion, setting collision factor to 0"<<endl;
//  log("samplingStatus")<<"Sampling ...\n";
//
//
//  double bias = 1.0; /// first round downhill, later use bias options.biasToTarget; /// bias to target between 0 and 1 (full bias)
//
//  gsl_vector* gradient = gsl_vector_alloc(protein->m_spanningTree->getNumDOFs());
//
//  std::list<Configuration*> samples;
//  samples.push_back(new Configuration(protein));
//
//  for(int i=0;i<options.samplesToGenerate;i++){
//    double dist = dist_to_objective(goal_distances);
//    if (dist<options.convergeDistance) break;
//
//    Configuration* seed = samples.back();
//    direction->gradient(seed, nullptr, gradient); //directed move
//
//    scale_gradient(gradient, options.maxRotation);
//    Configuration* new_conf = move->move(seed, gradient);
//    if (new_conf->updatedMolecule()->inCollision()) {
//      cout<<"Could not find new conformation"<<endl;
//      delete new_conf;
//      new_conf = nullptr;
//    }
//    else {
//      samples.push_back(new_conf);
//      log("samplingStatus") << "> New structure: conf_" + std::to_string((long long) i) + ".pdb" << endl;
//      log("samplingStatus") << "Distance to objective: " << dist_to_objective(goal_distances) << endl;
//      string fname = "output/conf_" + std::to_string((long long) i) + ".pdb";
//      IO::writePdb(new_conf->updatedMolecule(), fname);
//    }
//  }
//  double end_time = timer.ElapsedTime();
//
//  double finalDist = dist_to_objective(goal_distances);
//  if ( finalDist < options.convergeDistance){
//    log("samplingStatus") << "Reached DEER distances down to "<<finalDist<<" in "<<end_time - start_time<<" seconds. "<<endl;
//  }
//  timer.Reset();
//  start_time = timer.LastElapsedTime();
//  ////%%%%%%%%%%%%%%%%%%%%%%DOWNHILL TO TARGET%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ///Option with using a "real" target protein
  //Initialize target on the final conformation
//  Selection resNetworkTarget(options.residueNetwork);
//  Molecule* target = protein->deepClone();
//
//  ///Now we have all constraints and desired roots present to build the tree
//  target->initializeTree(resNetworkTarget,options.collisionFactor,options.roots);
//
//  if(options.collapseRigid>0) {
//    target = target->collapseRigidBonds(options.collapseRigid);
//  }
//
//  // Reset the protein
//  protein->setConfiguration(samples.front());
 */

  // Initialize planner
  DEERPlanner* planner;
  if(options.planner_string=="deer")
    planner = new DEERPlanner(
        protein,
        randDirection,
        move,
        direction,
        goal_distances,
        options.collisionCheck,
        blendedDir,
        options.samplesToGenerate,
        options.frontSize,
        options.stepSize,
        options.convergeDistance,
        options.biasToTarget
    );
  else{
    cerr<<"Unsuitable planner option for this application specified!"<<endl;
    exit(-1);
  }
  planner->initialize(randMove, metric, options.workingDirectory, options.saveData);

  //Use the planner "wrapper" function generateSamples() with its 3 stages approach, navigate, explore
  planner->generateSamples();

  //Print final status
  double end_time = timer.ElapsedTime();

  std::list<Configuration*>& m_samples = planner->getSamples();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianAndNullspaceTime<<" seconds\n";
  log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
  log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";

  if(options.saveData > 0){
    log("samplingStatus")<<"Creating trajectory"<<endl;
  }
  planner->createTrajectory();

  log("samplingStatus")<<"Done with approaching, start exploring"<<endl<<endl;

  ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New planner, protein with DEER distances, explore around?
  int remainingSamples = options.samplesToGenerate - (m_samples.size()-1);
  if(options.explore && remainingSamples > 0){
    Configuration* bestConfig = planner->getClosestConfig();
    protein->setConfiguration(bestConfig);

    Molecule* exploreProtein = protein->deepClone();
    const std::string name = "explored";
    exploreProtein->setName(name);

    ///Fixing DEER distances
    for (auto &goalDistance : goal_distances) {
      Atom *a1 = exploreProtein->getAtom(get<0>(goalDistance)->getId());
      Atom *a2 = exploreProtein->getAtom(get<1>(goalDistance)->getId());
      log("samplingStatus")<< "Adding distance constraint between "<<a1->getId()<<" and "<<a2->getId()<<endl;
      DBond *new_db = new DBond(a1, a2);
      exploreProtein->addDBond(new_db);
    }

    exploreProtein->initializeTree(resNetwork,options.collisionFactor,options.roots);

    //Clean-up old stuff to make room
    delete protein;
    delete planner;
    delete move;
    delete direction;

    if(options.collapseRigid>0) {
      exploreProtein = exploreProtein->collapseRigidBonds(options.collapseRigid);
      exploreProtein->writeRigidbodyIDToBFactor();
      string fname = "output/collapsed_explore.pdb";
      string fnamePml = "output/collapsed_explore.pml";
      IO::writePdb(exploreProtein, fname);
      IO::writePyMolScript(exploreProtein,fname,fnamePml);
    }

    log("samplingStatus")<<"DEER restrained collapsed molecule has:"<<endl;
    log("samplingStatus")<<"> "<<exploreProtein->m_spanningTree->m_cycleAnchorEdges.size()<<" total bond constraints"<<endl;
    log("samplingStatus")<<"> "<<exploreProtein->getDBonds().size()<<" distance bonds"<<endl;
    log("samplingStatus")<<"> "<<exploreProtein->m_spanningTree->getNumDOFs() << " DOFs of which " << exploreProtein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;


    auto explorePlanner = new DihedralRRT(exploreProtein,
                                          randDirection,
                                          remainingSamples,
                                          options.explorationRadius,
                                          options.maxRotation,
                                          false);

    explorePlanner->initialize(randMove, metric, options.workingDirectory, options.saveData);

    timer.Reset();
    start_time = timer.LastElapsedTime();

    explorePlanner->generateSamples();

    end_time = timer.ElapsedTime();

    std::list<Configuration*>& exploreSamples = explorePlanner->getSamples();
    log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(exploreSamples.size()-1)<<" valid samples\n";

    if(options.saveData > 0){
      log("samplingStatus")<<"Creating trajectory to most distant structure"<<endl;
    }
    explorePlanner->createTrajectory();

    delete explorePlanner;
    delete exploreProtein;
    ///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  }
  else {

    //Clean up
    delete protein;
    delete planner;
    delete move;
    delete direction;
  }

  delete randDirection;
  delete randMove;

  reportStream.close();
  plannerStream.close();
  debugStream.close();

  return 0;
}


void scale_gradient(gsl_vector* gradient, double maxRotation)
{
  double max=1;
  for(int i=0;i<gradient->size;i++){
    double val = gsl_vector_get(gradient, i);
    //double maxval = dof->getMaxPerturbation();
    if(fabs(val/maxRotation)>max)
      max = fabs(val/maxRotation);
  }
  gsl_vector_scale(gradient, 1/max);
}

double dist_to_objective(std::vector< std::tuple<Atom*, Atom*, double> > &goal_distances)
{
  double d=0;
  for (int i=0; i< goal_distances.size();i++){
    Coordinate c1= get<0>(goal_distances[i])->m_position;
    Coordinate c2= get<1>(goal_distances[i])->m_position;

//    d=d+fabs((sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z))-get<2>(goal_distances[i])));
    double dterm = sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z))-get<2>(goal_distances[i]);
    d += dterm*dterm;
  }
//  return d;
  return sqrt(d);
}


void printStrain(const Molecule& mol, const DeerOptions& options)
{

}