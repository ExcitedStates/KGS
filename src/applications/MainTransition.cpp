#include <stdexcept>
#include <string>
#include <iostream>
#include <list>
#include <math/NullspaceSVD.h>
#include <moves/ClashAvoidingMove.h>
#include <planners/SamplingPlanner.h>
#include <planners/PoissonPlanner.h>
#include <planners/BidirectionalMovingFront.h>
#include <planners/RRTPlanner.h>
#include <planners/DihedralRRT.h>

#include "core/Grid.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "moves/NullspaceMove.h"
#include "directions/RandomDirection.h"
#include "directions/DihedralDirection.h"
#include "directions/MSDDirection.h"
#include "directions/LSNullspaceDirection.h"
#include "directions/BlendedDirection.h"
#include "moves/DecreaseStepMove.h"
#include "metrics/RMSDnosuper.h"
#include "applications/options/TransitionOptions.h"

extern double jacobianAndNullspaceTime;
extern double rigidityTime;
extern double selectNodeTime;

using namespace std;

void targetedSampling(TransitionOptions& options){

  string pdb_file = options.initialStructureFile;
  Selection resNetwork(options.residueNetwork);
  Molecule* protein = IO::readPdb( pdb_file, resNetwork, options.extraCovBonds );
  protein->setCollisionFactor(options.collisionFactor);
  log() << "Molecule has " << protein->getAtoms().size() << " atoms\n";

  if(!options.annotationFile.empty())
    IO::readAnnotations(protein, options.annotationFile);

  // Do the same for the target
  string target_file = options.targetStructureFile;
  Molecule* target = IO::readPdb( target_file, resNetwork, options.extraCovBonds );
  target->setCollisionFactor(options.collisionFactor);

  //makes sure we have the same hydrogen bonds in target and m_molecule (m_molecule hbonds is adapted as well)
//  target->setToHbondIntersection(protein); // if desired, this has to be moved to before the tree construction

  /// Rigid bodies, spanning trees, and initial collisions
//  options.setResidueNetwork(protein);
//  options.setAtomSets(protein,target);

//  IO::readRigidbody( protein, resNetwork );
//  protein->buildRigidBodies();

//  unsigned int bestProteinRBId = protein->findBestRigidBodyMatch(options.m_root);//Todo: adapt this to usage without target
//  protein->buildSpanningTree(bestProteinRBId, options.flexibleRibose);//with the rigid body tree in place, we can generate a configuration
  //TODO: With multi-chain the choice of chain roots must be redesigned or removed
//  protein->buildSpanningTree(options.roots);//with the rigid body tree in place, we can generate a configuration

  /// Done in IO right now; ToDo: write wrapper, move out of IO
//  protein->setConfiguration(new Configuration(protein));

  // Check for collision
  // This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
//  protein->m_initialCollisions = protein->getAllCollisions();
  for(auto const& coll: protein->getInitialCollisions()){
    log("dominik")<<"Ini coll: "<<coll.first->getId()<<" "<<coll.first->getName()<<" "<<coll.second->getId()<<coll.second->getName()<<endl;
  }

  /// Rigid bodies, spanning trees, and initial collisions for the target
//  IO::readRigidbody( target, resNetwork );
//  target->buildRigidBodies();
  //Build rigid body tree for target
//  unsigned int bestTargetRBId = target->findBestRigidBodyMatch(options.m_root, &protein);
//  target->buildSpanningTree(bestTargetRBId, options.flexibleRibose);
//  target->buildSpanningTree(options.roots);


//  target->setConfiguration(new Configuration(target));

  // Check for collision
//  target->m_initialCollisions = target->getAllCollisions();
//    	for(mit=target->m_initialCollisions.begin(); mit != target->m_initialCollisions.end();mit++){
//    		Atom* atom1=mit->second.first;
//        	log("dominik")<<"Ini coll target: "<< mit->first.first << " "<< mit->second.first->getName() << " " << mit->first.second << mit->second.second->getName() <<endl;
//    	}


//	m_molecule.m_spanningTree->print();
  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus") << "> " << protein->getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus") << "> " << protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  log("samplingStatus")<<"Target has:"<<endl;
  log("samplingStatus")<<"> "<<target->getAtoms().size()<<" atoms"<<endl;
  log("samplingStatus")<<"> "<<target->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<target->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<target->m_spanningTree->getNumDOFs()<<" DOFs of which "<<target->m_spanningTree->getNumCycleDOFs()<<" are cycle-DOFs\n"<<endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  Selection metricSelection(options.metricSelection);
  try {
    if(TransitionOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(TransitionOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(TransitionOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Alignment
  if(options.alignIni){
    Selection alignSelection(options.alignSelection);
    //Alignment invalidates the current configuration, so reset new configuration
    Configuration* iniConf = target->m_conf;
    double initialRMSD = target->alignReferencePositionsTo(protein,alignSelection);//backup the aligned configuration
    target->setConfiguration(iniConf);
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
    move = new NullspaceMove(TransitionOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  //Initialize m_direction
  Direction* direction;
  bool blendedDir = false;
  Selection blendedSelection("all");
  Selection gradientSelection(options.gradientSelection);
  if(options.gradient == 0)
    direction = new RandomDirection(gradientSelection);
  else if(options.gradient == 1)
    direction = new DihedralDirection(gradientSelection);
  else if(options.gradient == 2){
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(new DihedralDirection(gradientSelection),0);
    bdir->addDirection(new RandomDirection(blendedSelection,TransitionOptions::getOptions()->maxRotation), 1);
    direction = bdir;
    blendedDir = true;
  }
  else if(options.gradient == 3)
    direction = new MSDDirection(gradientSelection, TransitionOptions::getOptions()->alignAlways);
  else if(options.gradient == 4){
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(new MSDDirection(gradientSelection, TransitionOptions::getOptions()->alignAlways),0);
    bdir->addDirection(new RandomDirection(blendedSelection,TransitionOptions::getOptions()->maxRotation), 1);
    direction = bdir;
    blendedDir = true;
  }
  else if(options.gradient <= 5)
    direction = new LSNullspaceDirection(gradientSelection);

  //Initialize planner
  SamplingPlanner* planner;
  if(options.planner_string=="binnedrrt")
    planner = new RRTPlanner(
        protein,
        direction,
        options.explorationRadius,
        options.samplesToGenerate,
        options.gradient,
        options.stepSize,
        options.maxRotation,
        options.scaleToRadius
    );
  else if(options.planner_string=="dihedralrrt")
    planner = new DihedralRRT(
        protein,
        direction,
        options.samplesToGenerate,
        options.explorationRadius,
        options.maxRotation,
        options.sampleRandom
    );
  else if(options.planner_string=="poisson")
    planner = new PoissonPlanner(
        protein,
        options.samplesToGenerate,
        options.poissonMaxRejectsBeforeClose,
        options.stepSize,
        options.gradientSelection
    );
  else if(options.planner_string=="dccrrt")
    planner = new BidirectionalMovingFront(
        protein,
        direction,
        target,
        metricSelection,
        options.collisionCheck,
        blendedDir,
        options.samplesToGenerate,
        options.frontSize,
        options.stepSize,
        options.switchAfter,
        options.convergeDistance,
        options.alignAlways,
        options.biasToTarget
    );
  else{
    cerr<<"Unknown planner option specified!"<<endl;
    exit(-1);
  }
  planner->initialize(move, metric, options.workingDirectory, options.saveData);

  if(options.saveData > 0){

    std::string out = options.workingDirectory + "output/" + target->getName() + "_lengths";
    IO::writeBondLengthsAndAngles(target, out);
    out = options.workingDirectory + "output/" + protein->getName() + "_lengths";
    IO::writeBondLengthsAndAngles(protein, out);

    if(options.saveData > 1){
      out = options.workingDirectory + "output/" + protein->getName() + "_q_target.txt";
      IO::writeQ(target,target->m_conf, out);
      out = options.workingDirectory + "output/" + protein->getName() + "_q_iniTarget.txt";
      IO::writeQ(target,protein->m_conf, out);
    }

    log() << "Total DOFs: " << protein->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs()
          << ", Max accessible DOFs: " << protein->m_spanningTree->getNumDOFs() - protein->m_spanningTree->getNumCycleDOFs() +
                                          protein->m_conf->getNullspace()->getNullspaceSize()
          << " Nullspace DOFs: " << protein->m_conf->getNullspace()->getNullspaceSize() << endl;fflush(stdout);
    log() << "Total DOFs in target: " << target->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << target->m_spanningTree->getNumCycleDOFs()
          << ", Max accessible DOFs: " << target->m_spanningTree->getNumDOFs() - target->m_spanningTree->getNumCycleDOFs() +
                                          target->m_conf->getNullspace()->getNullspaceSize()
          << " Nullspace DOFs: " << target->m_conf->getNullspace()->getNullspaceSize() << endl;fflush(stdout);

    if(options.saveData > 1){
      string out = options.workingDirectory + "output/" + protein->getName() + "_q_0.txt";
      IO::writeQ(protein, protein->m_conf, out);
    }

    log()<<"Number of rigid clusters: "<<protein->m_conf->m_numClusters;
    log()<<", biggest cluster: index "<<protein->m_conf->m_maxIndex<<" with "<<protein->m_conf->m_maxSize<<" atoms!"<<endl;
    //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigid << " rigidified and " << m_molecule.m_conf->CycleNullSpace->m_numCoordinated << " coordinated dihedrals" <<endl;
    //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigidHBonds<<" rigid out of "<<m_molecule.H_bonds.size()<<" hydrogen bonds!"<<endl<<endl;
    log()<< protein->m_conf->getNullspace()->getNumRigidDihedrals() << " rigidified";
    log()<<" and " << ( protein->m_conf->getNullspace()->getNumDOFs()-
                        protein->m_conf->getNullspace()->getNumRigidDihedrals()) << " coordinated dihedrals" <<endl;
    log()<< protein->m_conf->getNullspace()->getNumRigidHBonds()<<" rigid out of "<<protein->getHBonds().size()<<" hydrogen bonds!"<<endl;


    log()<<"Initial Distance: "<<metric->distance(protein->m_conf,target->m_conf)<<endl;

    log("samplingStatus")<<"Sampling ...\n"<<endl;
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
//    log("samplingStatus")<< planner->initialRebuildsAccepted<<" initial rebuild perturbations accepted (RI_ACC)"<<endl;
//    log("samplingStatus")<< planner->initialRebuildsRejected<<" initial rebuild perturbations rejected for clashing (RI_CLASH)"<<endl;
//    log("samplingStatus")<< planner->rebuildsAccepted<<" rebuild perturbations accepted (R_ACC)"<<endl;
//    log("samplingStatus")<< planner->rebuildsRejected<<" rebuild perturbations rejected for clashing (R_CLASH)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sAccepted<<" nullspace perturbations accepted (NS_ACC)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sRejected<<" nullspace perturbations rejected for clashing (NS_CLASH)"<<endl;


    if(options.saveData > 0){
      log("samplingStatus")<<"Creating trajectory"<<endl;
    }
    planner->createTrajectory();
  }
  log("samplingStatus")<<"Done"<<endl;

  //Clean up
  delete planner;
  delete move;
  delete direction;
  delete target;
  delete protein;
}


int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

//  ofstream reportStream;
//  reportStream.open("kgs_report.log");
//  enableLogger("report", reportStream);

//  ofstream plannerStream;
//  plannerStream.open("kgs_planner.log");
//  enableLogger("dominik", plannerStream);

//  ofstream debugStream;
//  debugStream.open("kgs_debug.log");
//  enableLogger("debug", debugStream);

  TransitionOptions::createOptions(argc, argv);

  TransitionOptions &options = *(TransitionOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//TransitionOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  // Set SVD cutoff
  NullspaceSVD::setSingularValueTolerance(options.svdCutoff);

  // Do the same for the target
  string target_file = options.targetStructureFile;
  if (target_file.empty()) {
    cerr << "MainTransition.cpp:" << __LINE__ << " ERROR: Must supply a target structure" << endl;
    return -1;
  }

  targetedSampling(options);

//  reportStream.close();
//  plannerStream.close();

  return 0;
}
