
#include <string>
#include <iostream>
#include <list>
#include <stdexcept>
#include "core/Molecule.h"
#include "core/Chain.h"
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
#include <moves/RawMove.h>
#include <moves/NullspaceMove.h>
#include <moves/ClashAvoidingMove.h>
#include "moves/CompositeMove.h"
#include <metrics/Metric.h>
#include <metrics/Dihedral.h>
#include <metrics/RMSD.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <directions/BlendedDirection.h>
#include <moves/DecreaseStepMove.h>
#include <metrics/RMSDnosuper.h>
#include <planners/PoissonPlanner2.h>

using namespace std;

extern double jacobianTime;
extern double rigidityTime;
extern double selectNodeTime;

void randomSampling(SamplingOptions& options); ///< randomized exploration without target structure
void targetedSampling(SamplingOptions& options); ///< directed/randomized exploration with target structure

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

  SamplingOptions::createOptions(argc, argv);

  SamplingOptions &options = *(SamplingOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//SamplingOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  // Do the same for the target
  string target_file = options.targetStructureFile;
  if (target_file.empty()) {
    cout<<"Random sampling without target"<<endl;
    randomSampling(options);
  }
  else {
    cout<<"Targeted sampling"<<endl;
    targetedSampling(options);
  }

  reportStream.close();
  plannerStream.close();

  return 0;
}

void randomSampling(SamplingOptions& options){

  string pdb_file = options.initialStructureFile;
  Molecule protein;
  protein.setCollisionFactor(options.collisionFactor);

  IO::readPdb( &protein, pdb_file, options.extraCovBonds );
  log() << "Molecule has " << protein.getAtoms().size() << " atoms\n";
  //options.setResidueNetwork(&protein);

  if(!options.annotationFile.empty())
    IO::readAnnotations(&protein, options.annotationFile);

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

  Selection resNetwork(options.residueNetwork);
  IO::readRigidbody( &protein, resNetwork );
//  unsigned int bestProteinRBId = protein.findBestRigidBodyMatch(options.m_root);//Todo: adapt this to usage without target
//  protein.buildSpanningTree(bestProteinRBId, options.flexibleRibose);//with the rigid body tree in place, we can generate a configuration
  //TODO: With multi-chain the choice of chain roots must be redesigned or removed
  protein.buildSpanningTree();//with the rigid body tree in place, we can generate a configuration
  cout<<"MainKGS .. spanning tree:"<<endl;

  protein.setConfiguration(new Configuration(&protein));

  // Check for collision
  // This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
  protein.m_initialCollisions = protein.getAllCollisions();
  for(auto const& coll: protein.m_initialCollisions){
    log("dominik")<<"Ini coll: "<<coll.first->getId()<<" "<<coll.first->getName()<<" "<<coll.second->getId()<<coll.second->getName()<<endl;
  }

//	m_molecule.m_spanning_tree->print();
  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein.getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein.m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->getNumDOFs() << " DOFs of which " << protein.m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  Selection metricSelection(options.metricSelection);
  try {
    if(SamplingOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new ClashAvoidingMove();
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(SamplingOptions::getOptions()->maxRotation);

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
      direction = new MSDDirection(selectionMoving, SamplingOptions::getOptions()->alignAlways);
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
    planner = new RRTPlanner(    &protein, *move, *metric, *direction );
  }else if(options.planner_string=="dihedralrrt"){
    log("samplingStatus")<<"Using dihedral RRT planner"<<endl;
    planner = new DihedralRRT(   &protein, *move, *metric, *direction );
  }else if(options.planner_string=="poisson"){
    log("samplingStatus")<<"Using Poisson-disk planner"<<endl;
    planner = new PoissonPlanner(&protein, *move, *metric );
  }else{
    cerr<<"Unknown planner option specified!"<<endl;
    exit(-1);
  }

  log() << "Total DOFs: " << protein.m_spanning_tree->getNumDOFs() << ", Cycle DOFs: " << protein.m_spanning_tree->getNumCycleDOFs()
        << ", Max accessible DOFs: " << protein.m_spanning_tree->getNumDOFs() - protein.m_spanning_tree->getNumCycleDOFs() + protein.m_conf->getNullspace()->NullspaceSize() << endl;fflush(stdout);

  if(options.saveData > 1){
    string out = options.workingDirectory + "output/" + protein.getName() + "_q_0.txt";
    IO::writeQ(&protein, protein.m_conf, out);
  }

  log()<<"Number of rigid clusters: "<<protein.m_conf->m_numClusters;
  log()<<", biggest cluster: index "<<protein.m_conf->m_maxIndex<<" with "<<protein.m_conf->m_maxSize<<" atoms!"<<endl;
  //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigid << " rigidified and " << m_molecule.m_conf->CycleNullSpace->m_numCoordinated << " coordinated dihedrals" <<endl;
  //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigidHBonds<<" rigid out of "<<m_molecule.H_bonds.size()<<" hydrogen bonds!"<<endl<<endl;
  log()<<protein.m_conf->getNullspace()->NumRigidDihedrals() << " rigidified";
  log()<<" and " << ( protein.m_conf->getNullspace()->getNumDOFs()-protein.m_conf->getNullspace()->NumRigidDihedrals()) << " coordinated dihedrals" <<endl;
  log()<<protein.m_conf->getNullspace()->NumRigidHBonds()<<" rigid out of "<<protein.getHBonds().size()<<" hydrogen bonds!"<<endl<<endl;

  log("samplingStatus")<<"Sampling ...\n";
  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  //Start exploring
  planner->GenerateSamples();

  //Print final status
  double end_time = timer.ElapsedTime();
  std::list<Configuration*> m_samples = planner->Samples();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianTime<<" seconds\n";
  log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
  log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";
//    log("samplingStatus")<< planner->initialRebuildsAccepted<<" initial rebuild perturbations accepted (RI_ACC)"<<endl;
//    log("samplingStatus")<< planner->initialRebuildsRejected<<" initial rebuild perturbations rejected for clashing (RI_CLASH)"<<endl;
//    log("samplingStatus")<< planner->rebuildsAccepted<<" rebuild perturbations accepted (R_ACC)"<<endl;
//    log("samplingStatus")<< planner->rebuildsRejected<<" rebuild perturbations rejected for clashing (R_CLASH)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sAccepted<<" nullspace perturbations accepted (NS_ACC)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sRejected<<" nullspace perturbations rejected for clashing (NS_CLASH)"<<endl;


  if(options.saveData > 0){
    planner->createTrajectory();
  }

  //Clean up
  delete planner;
  delete direction;
}

void targetedSampling(SamplingOptions& options){

  string pdb_file = options.initialStructureFile;
  Molecule protein;
  protein.setCollisionFactor(options.collisionFactor);

  IO::readPdb( &protein, pdb_file, options.extraCovBonds );
  log() << "Molecule has " << protein.getAtoms().size() << " atoms\n";

  if(!options.annotationFile.empty())
    IO::readAnnotations(&protein, options.annotationFile);

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

  // Do the same for the target
  string target_file = options.targetStructureFile;
  Molecule * target = new Molecule();
  target->setCollisionFactor(options.collisionFactor);
  IO::readPdb( target, target_file, options.extraCovBonds,&protein);

  //makes sure we have the same hydrogen bonds in target and m_molecule (m_molecule hbonds is adapted as well)
  target->setToHbondIntersection(&protein);

  /// Rigid bodies, spanning trees, and initial collisions
//  options.setResidueNetwork(&protein);
//  options.setAtomSets(&protein,target);

  Selection resNetwork(options.residueNetwork);
  IO::readRigidbody( &protein, resNetwork );

//  unsigned int bestProteinRBId = protein.findBestRigidBodyMatch(options.m_root);//Todo: adapt this to usage without target
//  protein.buildSpanningTree(bestProteinRBId, options.flexibleRibose);//with the rigid body tree in place, we can generate a configuration
  //TODO: With multi-chain the choice of chain roots must be redesigned or removed
  protein.buildSpanningTree();//with the rigid body tree in place, we can generate a configuration

  protein.setConfiguration(new Configuration(&protein));

  // Check for collision
  // This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
  protein.m_initialCollisions = protein.getAllCollisions();
  for(auto const& coll: protein.m_initialCollisions){
    log("dominik")<<"Ini coll: "<<coll.first->getId()<<" "<<coll.first->getName()<<" "<<coll.second->getId()<<coll.second->getName()<<endl;
  }

  /// Rigid bodies, spanning trees, and initial collisions for the target
  IO::readRigidbody( target, resNetwork );
  //Build rigid body tree for target
//  unsigned int bestTargetRBId = target->findBestRigidBodyMatch(options.m_root, &protein);
//  target->buildSpanningTree(bestTargetRBId, options.flexibleRibose);
  target->buildSpanningTree();

  //Alignment and spanning trees with possibly best m_root
  if(options.alignIni){
    target->alignReferencePositionsTo(&protein);//backup the aligned configuration
  }

  target->setConfiguration(new Configuration(target));

  // Check for collision
  target->m_initialCollisions = target->getAllCollisions();
//    	for(mit=target->m_initialCollisions.begin(); mit != target->m_initialCollisions.end();mit++){
//    		Atom* atom1=mit->second.first;
//        	log("dominik")<<"Ini coll target: "<< mit->first.first << " "<< mit->second.first->getName() << " " << mit->first.second << mit->second.second->getName() <<endl;
//    	}


//	m_molecule.m_spanning_tree->print();
  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus") << "> " << protein.getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein.m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus") << "> " << protein.m_spanning_tree->getNumDOFs() << " DOFs of which " << protein.m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  log("samplingStatus")<<"Target has:"<<endl;
  log("samplingStatus")<<"> "<<target->getAtoms().size()<<" atoms"<<endl;
  log("samplingStatus")<<"> "<<target->m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<target->m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<target->m_spanning_tree->getNumDOFs()<<" DOFs of which "<<target->m_spanning_tree->getNumCycleDOFs()<<" are cycle-DOFs\n"<<endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  Selection metricSelection(options.metricSelection);
  try {
    if(SamplingOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new ClashAvoidingMove();
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(SamplingOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  //Initialize direction
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
    bdir->addDirection(new RandomDirection(blendedSelection,SamplingOptions::getOptions()->maxRotation), 1);
    direction = bdir;
    blendedDir = true;
  }
  else if(options.gradient == 3)
    direction = new MSDDirection(gradientSelection, SamplingOptions::getOptions()->alignAlways);
  else if(options.gradient == 4){
    BlendedDirection* bdir = new BlendedDirection();
    bdir->addDirection(new MSDDirection(gradientSelection, SamplingOptions::getOptions()->alignAlways),0);
    bdir->addDirection(new RandomDirection(blendedSelection,SamplingOptions::getOptions()->maxRotation), 1);
    direction = bdir;
    blendedDir = true;
  }
  else if(options.gradient <= 5)
    direction = new LSNullspaceDirection(gradientSelection);

  //Initialize planner
  SamplingPlanner* planner;
  if(options.planner_string=="binnedrrt")         planner = new RRTPlanner(    &protein, *move, *metric, *direction );
  else if(options.planner_string=="dihedralrrt")  planner = new DihedralRRT(   &protein, *move, *metric, *direction );
  else if(options.planner_string=="poisson")      planner = new PoissonPlanner(&protein, *move, *metric);
  else if(options.planner_string=="dccrrt")       planner = new BidirectionalMovingFront(&protein, *move, *metric, *direction, target, blendedDir);
  else{
    cerr<<"Unknown planner option specified!"<<endl;
    exit(-1);
  }

  if(options.saveData > 0){
    std::string out = options.workingDirectory + "output/" + protein.getName() + "_target_lengths";
    IO::writeBondLengthsAndAngles(target, out);
    if(options.saveData > 1){
      out = options.workingDirectory + "output/" + protein.getName() + "_q_target.txt";
      IO::writeQ(target,target->m_conf, out);
      out = options.workingDirectory + "output/" + protein.getName() + "_q_iniTarget.txt";
      IO::writeQ(target,protein.m_conf, out);
    }

    log() << "Total DOFs: " << protein.m_spanning_tree->getNumDOFs() << ", Cycle DOFs: " << protein.m_spanning_tree->getNumCycleDOFs()
    << ", Max accessible DOFs: " << protein.m_spanning_tree->getNumDOFs() - protein.m_spanning_tree->getNumCycleDOFs() + protein.m_conf->getNullspace()->NullspaceSize() << endl;fflush(stdout);
    log() << "Total DOFs in target: " << target->m_spanning_tree->getNumDOFs() << ", Cycle DOFs: " << target->m_spanning_tree->getNumCycleDOFs()
    << ", Max accessible DOFs: " << target->m_spanning_tree->getNumDOFs() - target->m_spanning_tree->getNumCycleDOFs() + target->m_conf->getNullspace()->NullspaceSize() << endl;fflush(stdout);

    if(options.saveData > 1){
      string out = options.workingDirectory + "output/" + protein.getName() + "_q_0.txt";
      IO::writeQ(&protein, protein.m_conf, out);
    }

    log()<<"Number of rigid clusters: "<<protein.m_conf->m_numClusters;
    log()<<", biggest cluster: index "<<protein.m_conf->m_maxIndex<<" with "<<protein.m_conf->m_maxSize<<" atoms!"<<endl;
    //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigid << " rigidified and " << m_molecule.m_conf->CycleNullSpace->m_numCoordinated << " coordinated dihedrals" <<endl;
    //log()<<m_molecule.m_conf->CycleNullSpace->m_numRigidHBonds<<" rigid out of "<<m_molecule.H_bonds.size()<<" hydrogen bonds!"<<endl<<endl;
    log()<<protein.m_conf->getNullspace()->NumRigidDihedrals() << " rigidified";
    log()<<" and " << ( protein.m_conf->getNullspace()->getNumDOFs()-protein.m_conf->getNullspace()->NumRigidDihedrals()) << " coordinated dihedrals" <<endl;
    log()<<protein.m_conf->getNullspace()->NumRigidHBonds()<<" rigid out of "<<protein.getHBonds().size()<<" hydrogen bonds!"<<endl;


    log()<<"Initial Distance: "<<metric->distance(protein.m_conf,target->m_conf)<<endl;

    log("samplingStatus")<<"Sampling ...\n"<<endl;
    CTKTimer timer;
    timer.Reset();
    double start_time = timer.LastElapsedTime();

    //Start exploring
    planner->GenerateSamples();

    //Print final status
    double end_time = timer.ElapsedTime();
    std::list<Configuration*> m_samples = planner->Samples();
    log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
    log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianTime<<" seconds\n";
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
  delete target;
  delete direction;
}
