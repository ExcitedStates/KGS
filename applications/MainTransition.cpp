#include <stdexcept>
#include <string>
#include <iostream>
#include <list>
#include "core/Configuration.h"
#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"

#include <moves/NullspaceMove.h>
#include <moves/FastClashAvoidingMove.h>
#include <metrics/Dihedral.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <directions/BlendedDirection.h>
#include <moves/DecreaseStepMove.h>
#include <metrics/RMSDnosuper.h>

using namespace std;

extern double jacobianTime;
extern double rigidityTime;
extern double selectNodeTime;

void scale_gradient(gsl_vector* gradient, Molecule* mol);

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

  //SamplingOptions options(argc,argv);
  SamplingOptions::createOptions(argc, argv);

  SamplingOptions &options = *(SamplingOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//SamplingOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  string pdb_file = options.initialStructureFile;
  Molecule protein;
  protein.setCollisionFactor(options.collisionFactor);
  IO::readPdb( &protein, pdb_file, options.extraCovBonds );
//
//  if(!options.annotationFile.empty())
//    IO::readAnnotations(&protein, options.annotationFile);

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

  //Read the rigid body of the protein
  IO::readRigidbody( &protein );
  protein.buildSpanningTree();

  protein.setConfiguration(new Configuration(&protein));

  protein.m_initialCollisions = protein.getAllCollisions();

  string target_pdb_file = options.targetStructureFile;
  Molecule target;
  target.setCollisionFactor(options.collisionFactor);
  IO::readPdb( &target, target_pdb_file, options.extraCovBonds, &protein);

//  if(options.hydrogenbondMethod=="user")
//    IO::readHbonds( &target, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="rnaview")
//    IO::readHbonds_rnaview( &target, options.hydrogenbondFile, options.annotationFile.empty() );
//  else if(options.hydrogenbondMethod=="first" || options.hydrogenbondMethod=="FIRST")
//    IO::readHbonds_first( &target, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="kinari" || options.hydrogenbondMethod=="KINARI")
//    IO::readHbonds_kinari( &target, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="vadar")
//    IO::readHbonds_vadar( &target, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="dssr")
//    IO::readHbonds_dssr( &target, options.hydrogenbondFile );

  //makes sure we have the same hydrogen bonds in target and protein (protein hbonds is adapted as well)
  target.setToHbondIntersection(&protein);

  //Read the rigid body of the protein
  IO::readRigidbody( &target );
  target.buildSpanningTree();

  target.setConfiguration(new Configuration(&target));

  target.m_initialCollisions = target.getAllCollisions();

//	m_molecule.m_spanning_tree->print();
  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein.getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein.m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->getNumDOFs() << " DOFs of which " << protein.m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  //Build rigid body tree for target
//  unsigned int bestTargetRBId = target->findBestRigidBodyMatch(options.m_root, &protein);
//  target->buildSpanningTree(bestTargetRBId, options.flexibleRibose);
  log("samplingStatus")<<"Target has:"<<endl;
  log("samplingStatus")<<"> "<<target.getAtoms().size()<<" atoms"<<endl;
  log("samplingStatus")<<"> "<<target.m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<target.m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<target.m_spanning_tree->getNumDOFs()<<" DOFs of which "<<target.m_spanning_tree->getNumCycleDOFs()<<" are cycle-DOFs\n"<<endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  try {
    Selection metricSelection(options.metricSelection);
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
    move = new FastClashAvoidingMove();
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(SamplingOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  Selection resNetwork(options.residueNetwork);

  //Initialize direction
  Direction* direction;
  if(options.gradient == 0)
    direction = new RandomDirection(resNetwork);
  else if(options.gradient == 1)
    direction = new DihedralDirection(resNetwork);
  else if(options.gradient == 2){
    BlendedDirection* m_direction = new BlendedDirection();
    m_direction->addDirection(new DihedralDirection(resNetwork),0);
    m_direction->addDirection(new RandomDirection(resNetwork,SamplingOptions::getOptions()->maxRotation), 1);
    direction = m_direction;
  }
  else if(options.gradient == 3)
    direction = new MSDDirection(resNetwork);
  else if(options.gradient == 4){
    BlendedDirection* m_direction = new BlendedDirection();
    m_direction->addDirection(new MSDDirection(resNetwork),0);
    m_direction->addDirection(new RandomDirection(resNetwork,SamplingOptions::getOptions()->maxRotation), 1);
    direction = m_direction;
  }else if(options.gradient == 5)
    direction = new LSNullspaceDirection(resNetwork);
  else{
    cerr<<__FILE__<<": Unknown gradient: "<<options.gradient<<endl;
    exit(-1);
  }


  if(options.saveData > 0){
    std::string out = options.workingDirectory + "output/" + protein.getName() + "_target_lengths";
    IO::writeBondLengthsAndAngles(&target, out);
    if(options.saveData > 1){
      out = options.workingDirectory + "output/" + protein.getName() + "_q_target.txt";
      IO::writeQ(&target,target.m_conf, out);
      out = options.workingDirectory + "output/" + protein.getName() + "_q_iniTarget.txt";
      IO::writeQ(&target,protein.m_conf, out);
    }

    log() << "Total DOFs: " << protein.m_spanning_tree->getNumDOFs() << ", Cycle DOFs: " << protein.m_spanning_tree->getNumCycleDOFs() << endl;fflush(stdout);
    log() << "Total DOFs in target: " << target.m_spanning_tree->getNumDOFs() << ", Cycle DOFs: " << target.m_spanning_tree->getNumCycleDOFs() << endl << endl;fflush(stdout);

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


    log()<<"Initial Distance: "<<metric->distance(protein.m_conf,target.m_conf);

    log("samplingStatus")<<"Sampling ...\n";
    CTKTimer timer;
    timer.Reset();
    double start_time = timer.LastElapsedTime();

    gsl_vector* gradient = gsl_vector_alloc(protein.m_spanning_tree->getNumDOFs());
    Configuration* target_conf = new Configuration(&target);
    std::list<Configuration*> samples;
    samples.push_back(new Configuration(&protein));
    for(int i=0;i<options.samplesToGenerate;i++){
      cout<<"Iteration "<<i<<endl;
      Configuration* seed = samples.back();
      direction->gradient(seed, target_conf, gradient);
      //gsl_vector_scale_max_component(gradient, options.maxRotation);
      scale_gradient(gradient, &protein);
      gsl_vector_scale(gradient, options.stepSize);
      Configuration* new_conf = move->move(seed, gradient);
      IO::writePdb(new_conf->updatedMolecule(), "output/conf_"+std::to_string((long long)i)+".pdb");
      samples.push_back(new_conf);
    }



    //Print final status
    double end_time = timer.ElapsedTime();
    log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(samples.size()-1)<<" valid samples\n";
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
  }
  log("samplingStatus")<<"Done"<<endl;
  //Clean up
  delete direction;


  return 0;
}


void scale_gradient(gsl_vector* gradient, Molecule* mol)
{
  double factor = 1.0;

  for(int i=0;i<mol->m_spanning_tree->getNumDOFs();i++){
    DOF* dof = mol->m_spanning_tree->getDOF(i);
    int idx = dof->getIndex();
    double val = gsl_vector_get(gradient, idx);
    double maxval = dof->getMaxPerturbation();
    if(fabs(maxval/val)<factor)
      factor = fabs(maxval/val);
  }
  gsl_vector_scale(gradient, factor);
}

