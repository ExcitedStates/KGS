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
#include <directions/BlendedDirection.h>
#include <math/gsl_helpers.h>
#include <gsl/gsl_vector_double.h>

using namespace std;

extern double selectNodeTime;

BidirectionalMovingFront::BidirectionalMovingFront(Molecule * protein, Move& move, metrics::Metric& metric, Direction& direction, Molecule * target,bool blendedDir):
    SamplingPlanner(move,metric),
    m_protein(protein),
    direction(direction),
    m_target(target),
    m_isBlended(blendedDir),
    m_rmsd(new metrics::RMSD(metric.getSelection())),
    stopAfter(SamplingOptions::getOptions()->samplesToGenerate),
    m_frontSize(SamplingOptions::getOptions()->frontSize)
{
  if (m_target == nullptr ){
    cerr<<"No valid target specified. Please provide target or choose a different planner."<<endl;
    exit(-1);
  }
  m_fwdRoot = new Configuration( m_protein );
  m_fwdRoot->updateMolecule();
  //m_fwdRoot->computeCycleJacobianAndNullSpace();
  m_fwdRoot->m_vdwEnergy = (m_protein->vdwEnergy(&(m_protein->m_initialCollisions),SamplingOptions::getOptions()->collisionCheck)).second;
  m_fwdSamples.push_back( m_fwdRoot );
  m_fwdFront.push_back( m_fwdRoot );

  m_revRoot = new Configuration( m_target );
  m_revRoot->updateMolecule();
  //m_revRoot->computeCycleJacobianAndNullSpace();
  m_revRoot->m_id = 1;
  m_revRoot->m_vdwEnergy = (m_target->vdwEnergy( &(m_target->m_initialCollisions),SamplingOptions::getOptions()->collisionCheck)).second;
  m_revSamples.push_back( m_revRoot );
  m_revFront.push_back( m_revRoot );

  //closest configs and distance measures
  m_minDistance = m_metric.distance(m_fwdRoot,m_revRoot);
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

std::list<Configuration*>& BidirectionalMovingFront::Samples(){

  std::list<Configuration*>& allSamples = m_fwdSamples;
  allSamples.merge(m_revSamples);
  return allSamples;
}

void BidirectionalMovingFront::GenerateSamples() {

  static int numSamples = 0, failedTrials = 0, totalTrials = 0;
  static double accept_ratio;
  double stepSize = SamplingOptions::getOptions()->stepSize;

  int samplesTillSwap = SamplingOptions::getOptions()->switchAfter;
  bool swapped = false;

  Configuration *qTarget = nullptr, *qSeed, *qNew = nullptr; //this qTarget is either global or random and can change for each sample

  while (numSamples < stopAfter) {//sample at most until the number of samples has been reached
    ++totalTrials;

    CTKTimer timer;
    double start_time = timer.getTimeNow();

    if (SamplingOptions::getOptions()->sampleRandom || qNew == nullptr ||
        m_addedToFront == false || swapped == true) {//create new target and compute closest seed
      log("dominik") << "Generating new target" << endl;
      swapped=false;
      if (qTarget != nullptr) {//prevent leakage, delete existing target configuration from previous round
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

    gsl_vector* gradient = gsl_vector_calloc(m_protein->totalDofNum());

    if(m_isBlended){
      BlendedDirection& blendedDir = reinterpret_cast<BlendedDirection&>(direction);
      blendedDir.changeWeight(0,double(numSamples)/double(stopAfter) );
      blendedDir.changeWeight(1,1.0 - double(numSamples)/double(stopAfter) );
    }
    direction.gradient(qSeed, qTarget, gradient); //computes the search direction for a new sample
    gsl_vector_scale_to_length(gradient, stepSize);
    qNew = move.move(qSeed, gradient); //Perform move

    if(qNew->updatedMolecule()->inCollision() ){
      failedTrials++;
      totalTrials++;
      delete qNew;
      qNew = nullptr;
      //Could not find a new conformation, swap search directions
      swapFwdRev();
    }
    else{//collision-free
      //Potentially reject new config if large violations?
      double violationNorm = m_protein->checkCycleClosure(qNew);

      numSamples++;
      samplesTillSwap--;
      gsl_vector_free(gradient);

      //Push-back in forward tree
      m_fwdSamples.push_back(qNew);
      qNew->m_id = numSamples;

      //Distance computations
      evaluateDistances(qNew);
      //Check if front needs to be updated
      updateFwdFront(qNew);

      //Enthalpy and entropy computation
//      log("dominik")<<"Length of all collisions: "<<allCollisions.size()<<endl;
//      pair<double,double> enthalpyVals=m_molecule->vdwEnergy(&allCollisions,m_options.collisionCheck);
//      qNew->m_vdw = enthalpyVals.first;
//      qNew->m_deltaH = enthalpyVals.second - new_q->m_vdw;

      qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(SamplingOptions::getOptions()->collisionCheck);

      log("dominik") << "> New structure: "<<qNew->getMolecule()->getName()<<"_new_"<<numSamples<<".pdb, accessible dofs: "<<qNew->m_clashFreeDofs<<endl<<endl;

      log("samplingStatus") << "> New structure: " << qNew->getMolecule()->getName()<<"_new_" << numSamples << ".pdb";
      log("samplingStatus") << " .. Distance to initial: " << setprecision(6) << qNew->m_distanceToIni;
      log("samplingStatus") << " .. Distance to moving-front target: " << setprecision(3) << qNew->m_paretoFrontDistance;
      log("samplingStatus") << " .. accessible dofs: "<<qNew->m_clashFreeDofs<<endl;
      log("samplingStatus") << " .. norm constraint violation: "<<violationNorm<<endl;
      log("samplingStatus")<<"Current shortest distance between both trees at configs "<<m_closestFwdSample->m_id <<" and "<<m_closestRevSample->m_id<<", Distance: "<<m_minDistance<<endl;

      writeNewSample(qNew, m_fwdRoot,qNew->m_id);

      //Check if target has been reached!
      if(m_minDistance < SamplingOptions::getOptions()->convergeDistance){
        log()<<"Reached target, not creating more samples!"<<" distance: "<<m_minDistance<<endl;
        break;
      }
      if(samplesTillSwap==0){
        swapFwdRev();
        samplesTillSwap = SamplingOptions::getOptions()->switchAfter;
        swapped = true;
      }
    }
  }

  log() << endl << "Max tree depth (excluding the m_root) = " << m_max_depth << endl;
  log() << "failed " << failedTrials << " times due to collision." << endl;

}

void BidirectionalMovingFront::updateFwdFront(Configuration *qNew) {

  std::list<Configuration*>::iterator cit;
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

  //Distance to target m_root
  m_revRoot->updateMolecule(); //go back to original, initial target and compute distance
  if(SamplingOptions::getOptions()->alignAlways){
    alignVal = m_rmsd->align(m_target, m_protein);
    qNew->m_distanceToTarget = m_rmsd->distance_noOptimization(qNew,m_revRoot);
  }else{
    qNew->m_distanceToTarget = m_metric.distance(qNew,m_revRoot);
  }

  //Distance to own m_root
  if(SamplingOptions::getOptions()->alignAlways) {
    qNew->m_distanceToIni = m_rmsd->distance_noOptimization(qNew, m_fwdRoot);
  }else{
    qNew->m_distanceToIni = m_metric.distance(qNew, m_fwdRoot);
  }

  //Distance to parents
  if(qNew->getParent()!= nullptr ) {
    qNew->m_distanceToParent = m_metric.distance(qNew, qNew->getParent());
  }

  m_closestRevSample->updateMolecule();//Distance to closest sample from opposing tree

  if(SamplingOptions::getOptions()->alignAlways){
    alignVal = m_rmsd->align(m_target, m_protein);
    qNew->m_paretoFrontDistance = m_rmsd->distance_noOptimization(qNew,m_closestRevSample);
  }else{
    qNew->m_paretoFrontDistance = m_metric.distance(qNew,m_closestRevSample);
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
    for (int i=0; i<pTarget->getNumDOFs(); ++i) {
      pTarget->m_dofs[i]=Math3D::dPi*RandomN1P1();
    }
    pTarget->m_id = -1; //invalid ID, identify non-sampled coonformations
  }

  pTarget->updateMolecule();

  return pTarget;
}

Configuration* BidirectionalMovingFront::SelectSeed (Configuration *pTarget) {

  Configuration *pSmp, *pMinSmp;
  double minDistance = 1000000.0;
  double distance;
  pMinSmp = nullptr;

//  const vector<int>* resNetwork = &(SamplingOptions::getOptions()->residueNetwork);
//  bool allResidues = resNetwork->size() == 0 ? true:false;
//  Selection resSelection(SamplingOptions::getOptions()->residueNetwork);
  bool allResidues = true; //TODO: Update so old code works again

  string rmsdSet="HEAVY", dihSet="MOV";
  if(!allResidues){
    rmsdSet="RESHEAVY";
    dihSet="RESMOV";
  }

  pTarget->updatedMolecule();//update target protein

  for (auto const& pSmp : m_fwdFront){
//    distance = m_metric.distance(pSmp,m_target->m_conf);//TODO: Why m_target->m_conf and not pTarget?
    distance = m_metric.distance(pSmp,pTarget);
    if (distance < minDistance) {
      minDistance = distance;
      pMinSmp = pSmp;
    }
  }

  return pMinSmp;
}

void BidirectionalMovingFront::swapFwdRev(){

  //Here, we set a sample from the existing forward front as target, and then switch search directions
  std::list<Configuration*>::iterator cit;
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
  //This will be the target for the next round
  log("dominik")<<"Using sample "<<pSmp->m_id<<" as current global target!"<<endl;
  pSmp->updateMolecule();
  m_currentGlobalTarget = pSmp;

  //Switch directions
  std::swap(m_fwdSamples,m_revSamples);
  std::swap(m_fwdFront,m_revFront);
  std::swap(m_fwdRoot,m_revRoot);
  std::swap(m_protein,m_target);
  std::swap(m_closestFwdSample,m_closestRevSample);

}


void BidirectionalMovingFront::createTrajectory(){

  m_closestFwdSample->updateMolecule();
  m_closestRevSample->updateMolecule();
  log()<<"The forward trajectory ends at sample "<<m_closestFwdSample->m_id<<endl;
  log()<<"The reverse trajectory ends at sample "<<m_closestRevSample->m_id<<endl;

  SamplingOptions& options = *(SamplingOptions::getOptions());

  const string& out_path = options.workingDirectory;
  const string& name = m_protein->getName();

  string out_collPdb = out_path + "output/" + name + "_path.pdb";
  ///save pyMol movie script
  string out_pyMol=out_path + "output/" +  name + "_pathMov.pml";

  IO::writeTrajectory(m_protein, out_collPdb, out_pyMol);
}

