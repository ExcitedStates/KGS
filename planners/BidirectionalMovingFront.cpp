//
// Created by Dominik Budday on 09.03.16.
//

#include "BidirectionalMovingFront.h"

#include "core/Molecule.h"
#include "core/Chain.h"
#include "Logger.h"
#include "core/Transformation.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "metrics/Dihedral.h"

#include <iomanip>
#include <stack>
#include <IO.h>

using namespace std;

extern double selectNodeTime;

BidirectionalMovingFront::BidirectionalMovingFront(Molecule * protein, Move& move, Direction& direction, Molecule * target):
    SamplingPlanner(move),
    m_protein(protein),
    direction(direction),
    m_target(target),
    stopAfter(SamplingOptions::getOptions()->samplesToGenerate),
    m_frontSize(SamplingOptions::getOptions()->frontSize)
{
  if (m_target == NULL ){
    cerr<<"No valid target specified. Please provide target or choose a different planner."<<endl;
    exit(-1);
  }
  m_fwdRoot = new Configuration( m_protein );
  m_fwdRoot->updateProtein();
  m_fwdRoot->computeCycleJacobianAndNullSpace();
  m_fwdRoot->m_vdwEnergy = (m_protein->vdwEnergy(&(m_protein->Initial_collisions),SamplingOptions::getOptions()->collisionCheck)).second;
  m_fwdSamples.push_back( m_fwdRoot );
  m_fwdFront.push_back( m_fwdRoot );

  m_revRoot = new Configuration( m_target );
  m_revRoot->updateProtein();
  m_revRoot->computeCycleJacobianAndNullSpace();
  m_revRoot->m_id = 1;
  m_revRoot->m_vdwEnergy = (m_target->vdwEnergy( &(m_target->Initial_collisions),SamplingOptions::getOptions()->collisionCheck)).second;
  m_revSamples.push_back( m_revRoot );
  m_revFront.push_back( m_revRoot );

  if(SamplingOptions::getOptions()->metric_string=="rmsd") 		  m_metric = new metrics::RMSD();
  if(SamplingOptions::getOptions()->metric_string=="dihedral") 	m_metric = new metrics::Dihedral();

  //closest configs and distance measures
  m_minDistance = m_metric->distance(m_fwdRoot,m_revRoot);
  m_closestFwdSample = m_fwdRoot;
  m_closestRevSample =  m_revRoot;
  m_currentGlobalTarget = m_revRoot;//target configuration for this iteration
  m_addedToFront = true;

  m_max_depth=0;

  m_nCDCall = 0;
  m_nMetricsCall = 1;

  movesRejected = 0;
  movesAccepted = 0;
}

BidirectionalMovingFront::~BidirectionalMovingFront() {
  Configuration *pSmp;
  for (auto iter=m_fwdSamples.begin(); iter!=m_fwdSamples.end(); iter++) {
    pSmp = *iter;
    delete pSmp;
  }
  for (auto iter=m_revSamples.begin(); iter!=m_revSamples.end(); iter++) {
    pSmp = *iter;
    delete pSmp;
  }
}

void BidirectionalMovingFront::GenerateSamples() {

  static int numSamples = 0, failedTrials = 0, totalTrials = 0;
  static double accept_ratio;

  Configuration *qTarget = NULL, *qSeed, *qNew = NULL; //this qTarget is either global or random and can change for each sample

  while (numSamples < stopAfter) {//sample at most until the number of samples has been reached
    ++totalTrials;

    CTKTimer timer;
    double start_time = timer.getTimeNow();

    if (SamplingOptions::getOptions()->sampleRandom || qNew == NULL ||
        m_addedToFront == false) {//create new target and compute closest seed
      log("dominik") << "Generating new target" << endl;

      if (qTarget != NULL) {//prevent leakage, delete existing target configuration from previous round
        delete qTarget;
      }

      qTarget = GenerateRandConf();//New target configuration: random or m_currentGlobalTarget, probability specified with bias
      qSeed = SelectSeed(qTarget); //Compute closest seed
      double end_time = timer.getTimeNow();
      selectNodeTime += end_time - start_time;
    }
    else {//we found a valid sample in the last round and aim towards the same target using the latest configuration as new seed
      log("dominik") << "Using latest sample and previous target" << endl;
      qSeed = m_fwdSamples.back();
    }

    log("dominik") << "Using sample " << qSeed->m_id << " as base" << endl;
    log("dominik") << "Using sample " << qTarget->m_id << " as target" << endl;

    gsl_vector* gradient = gsl_vector_alloc(m_protein->totalDofNum());
    direction.gradient(qSeed, qTarget, gradient); //computes the search direction for a new sample
    qNew = move.move(qSeed, gradient); //Perform move

    if(qNew->updatedProtein()->inCollision() ){
      failedTrials++;
      totalTrials++;
      delete qNew;
      qNew = NULL;
      //Could not find a new conformation, swap search directions
      swapFwdRev();
    }
    else{//collision-free
      //Potentially reject new config if large violations?
      m_protein->checkCycleClosure(qNew);

      numSamples++;
      gsl_vector_free(gradient);

      //Push-back in forward tree
      m_fwdSamples.push_back(qNew);
      qNew->m_id = numSamples;

      //Distance computations
      evaluateDistances(qNew);
      //Check if front needs to be updated
      updateFwdFront(qNew);

      //Enthalpy and entropy computation, currently in move
//      log("dominik")<<"Length of all collisions: "<<allCollisions.size()<<endl;
//      pair<double,double> enthalpyVals=m_protein->vdwEnergy(&allCollisions,m_options.collisionCheck);
//      new_q->m_vdw = enthalpyVals.first;
//      new_q->m_deltaH = enthalpyVals.second - new_q->m_vdw;

      log("samplingStatus") << "> New structure: "<<SamplingOptions::getOptions()->moleculeName<<"_new_"<<numSamples<<".pdb"<<endl;
      log()<<"Distance to target: "<<qNew->m_distanceToTarget<<", accessible dofs: "<<qNew->m_clashFreeDofs<<endl;
      log()<<"Current shortest distance between both trees at configs "<<m_closestFwdSample->m_id <<" and "<<m_closestRevSample->m_id<<", Distance: "<<m_minDistance<<endl;

      writeNewSample(qNew, m_fwdRoot,qNew->m_id);

      //Check if target has been reached!
      if(m_minDistance < SamplingOptions::getOptions()->convergeDistance){
        log()<<"Reached target, not creating more samples!"<<" distance: "<<m_minDistance<<endl;
        break;
      }
    }
  }

  log() << endl << "Max tree depth (excluding the root) = " << m_max_depth << endl;
  log() << "failed " << failedTrials << " times due to collision." << endl;

}

void BidirectionalMovingFront::updateFwdFront(Configuration *qNew) {

  ConfigurationList::iterator cit;
  int i=1;
  m_addedToFront = false;

  for( cit=m_fwdFront.begin(); cit!=m_fwdFront.end(); cit++ ){
    if(qNew->m_paretoFrontDistance <= (*cit)->m_paretoFrontDistance){
      log("dominik")<<"Inserting into pareto front at: "<<i<<", with distance: "<<qNew->m_paretoFrontDistance<<endl;
      m_fwdFront.insert(cit,qNew); //potentially use it as new seed / target for reverse planner
      m_addedToFront = true;
      break;
    }
    i++;
  }
  if(cit==m_fwdFront.end() && m_fwdFront.size() < m_frontSize){
    m_fwdFront.push_back(qNew);
    log("dominik")<<"Inserting at the end at: "<<i<<", with distance: "<<qNew->m_paretoFrontDistance<<endl;
    m_addedToFront = true;
  }
  if(m_fwdFront.size() > m_frontSize)
    m_fwdFront.pop_back();//keep length at maximum
}

void BidirectionalMovingFront::evaluateDistances(Configuration* qNew){

  double alignVal;

  //Distance to target root
  m_revRoot->updateProtein(); //go back to original, initial target and compute distance
  if(SamplingOptions::getOptions()->alignAlways){
    alignVal = metrics::RMSD::align(m_target, m_protein);
    qNew->m_distanceToTarget = metrics::RMSD::distance_noOptimization(qNew,m_revRoot);
  }else{
    qNew->m_distanceToTarget = m_metric->distance(qNew,m_revRoot);
  }

  //Distance to own root
  if(SamplingOptions::getOptions()->alignAlways) {
    qNew->m_distanceToIni = metrics::RMSD::distance_noOptimization(qNew, m_fwdRoot);
  }
  else{
    qNew->m_distanceToIni = m_metric->distance(qNew, m_fwdRoot);
  }
  //Distance to parents
  if(qNew->getParent()!= NULL && qNew->getParent()->getParent() != NULL) {
    double distToParent = m_metric->distance(qNew, qNew->getParent());
//    double distToParentsParent = m_metric->distance(qNew,qNew->getParent()->getParent());
    qNew->m_distanceToParent = distToParent;
  }

  m_closestRevSample->updateProtein();//Distance to closest sample from opposing tree

  if(SamplingOptions::getOptions()->alignAlways){
    alignVal = metrics::RMSD::align(m_target, m_protein);
    qNew->m_paretoFrontDistance = metrics::RMSD::distance_noOptimization(qNew,m_closestRevSample);
  }else{
    qNew->m_paretoFrontDistance = m_metric->distance(qNew,m_closestRevSample);

  }

  if( qNew->m_paretoFrontDistance < m_minDistance){//update current shortest configs
    m_minDistance = qNew->m_paretoFrontDistance;
    m_closestFwdSample = qNew;
  }

  if( qNew->m_treeDepth > m_max_depth)
    m_max_depth=qNew->m_treeDepth;
}

Configuration* BidirectionalMovingFront::GenerateRandConf() {

  double randVal = 1.0;
  //Compute the target configuration: either random or a sampled target conf (depending on bias)
  Configuration* pTarget;
  randVal = Random01();

  if(randVal <= SamplingOptions::getOptions()->biasToTarget ){//use directly the unperturbed target configuration
    log("dominik")<<"Using global target"<<endl;
    pTarget = m_currentGlobalTarget->clone();
  }
  else{
    pTarget = new Configuration(m_revRoot);
    log("dominik")<<"Using random target"<<endl;
    for (int i=0; i<pTarget->m_numDOFs; ++i) {
      pTarget->m_dofs[i]=dPi*RandomN1P1();
    }
    pTarget->m_id = -1; //invalid ID, identify non-sampled coonformations
  }

  pTarget->updateProtein();

  return pTarget;
}

Configuration* BidirectionalMovingFront::SelectSeed (Configuration *pTarget) {

  Configuration *pSmp, *pMinSmp;
  double minDistance = 1000000.0;
  double distance;
  pMinSmp = NULL;

  const vector<int>* resNetwork = &(SamplingOptions::getOptions()->residueNetwork);
  bool allResidues = resNetwork->size() == 0 ? true:false;

  string rmsdSet="HEAVY", dihSet="MOV";
  if(!allResidues){
    rmsdSet="RESHEAVY";
    dihSet="RESMOV";
  }

  pTarget->updatedProtein();//update target protein

  for (list<Configuration*>::iterator iter=m_fwdFront.begin(); iter!=m_fwdFront.end(); ++iter) {
    pSmp = *iter;
    if(SamplingOptions::getOptions()->alignAlways){
      double alignVal = metrics::RMSD::align(m_target, pSmp->updatedProtein());
      distance = m_metric->distance(pSmp,m_target->m_conf);
    }
    else{
      distance = m_metric->distance(pSmp,m_target->m_conf);
    }
    if (distance < minDistance) {
      minDistance = distance;
      pMinSmp = pSmp;
    }
  }

  return pMinSmp;
}

void BidirectionalMovingFront::swapFwdRev(){

  //Here, we set a sample from the existing forward front as target, and then switch search directions
  ConfigurationList::iterator cit;
  Configuration *pSmp;
  log("dominik")<<"Size of pareto front: "<<m_fwdFront.size()<<endl;
  int randConf = rand() % (m_fwdFront.size()); //todo: change this to more often use the "best sample"
  log("dominik")<<"Choosing rand conf number "<<randConf<<endl;
  int i=0;
  for( cit=m_fwdFront.begin(); cit!=m_fwdFront.end(); cit++ ){
    if(i == randConf){
      pSmp = (*cit);
      break;
    }
    i++;
  }
  if(cit == m_fwdFront.end() ){
    cerr<<"No valid target sample"<<endl;
    exit(-1);
  }
  //This will be the target for the next round
  log("dominik")<<"Using sample "<<pSmp->m_id<<" as current global target!"<<endl;
  pSmp->updateProtein();
  m_currentGlobalTarget = pSmp;

  //Switch directions
  std::swap(m_fwdSamples,m_revSamples);
  std::swap(m_fwdFront,m_revFront);
  std::swap(m_fwdRoot,m_revRoot);
  std::swap(m_protein,m_target);
  std::swap(m_closestFwdSample,m_closestRevSample);

}


void BidirectionalMovingFront::createTrajectory(){

  m_closestFwdSample->updateProtein();
  m_closestRevSample->updateProtein();
  log()<<"The forward trajectory ends at sample "<<m_closestFwdSample->m_id<<endl;
  log()<<"The reverse trajectory ends at sample "<<m_closestRevSample->m_id<<endl;

  SamplingOptions& options = *(SamplingOptions::getOptions());

  const string& out_path = options.workingDirectory;
  const string& name = options.moleculeName;

  string out_collPdb = out_path + "output/" + name + "_path.pdb";
  ///save pyMol movie script
  string out_pyMol=out_path + "output/" +  name + "_pathMov.pml";

  IO::writeTrajectory(m_protein, out_collPdb, out_pyMol);
}

