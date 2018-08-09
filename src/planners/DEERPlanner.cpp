#include <utility>

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
// Created by Dominik Budday on 31.07.18.
//

#include "DEERPlanner.h"

#include <iomanip>
#include <stack>

#include "Logger.h"
#include "core/Transformation.h"
#include "CTKTimer.h"
#include "metrics/Dihedral.h"

#include "IO.h"
#include "directions/BlendedDirection.h"
#include "math/gsl_helpers.h"
#include "core/DBond.h"

using namespace std;

extern double selectNodeTime;
//double dist_to_objective(std::vector< std::tuple<Atom*, Atom*, double> > &goal_distances);

DEERPlanner::DEERPlanner(
    Molecule *protein,
    Direction *direction,
    Move *directedMove,
    Direction *directedDirection,
    std::vector< std::tuple<Atom*, Atom*, double> > goal_distances,
    string& collisionCheck,
    bool blendedDir,
    int stopAfter,
    int frontSize,
    double stepSize,
    double convergeDistance,
    double biasToTarget
) :
    SamplingPlanner(),
    m_protein(protein),
    m_direction(direction),
    m_directedMove(directedMove),
    m_directedDirection(directedDirection),
    m_goalDistances(std::move(goal_distances)),
    m_isBlended(blendedDir),
    m_stopAfter(stopAfter),
    m_frontSize(frontSize),
    m_collisionCheck(collisionCheck),
    m_stepSize(stepSize),
    m_convergeDistance(convergeDistance),
    m_biasToTarget(biasToTarget)
{
  m_fwdRoot = m_protein->m_conf;
  m_fwdRoot->m_id = 0;
  m_fwdRoot->m_vdwEnergy = (m_protein->vdwEnergy(&(m_protein->getInitialCollisions()), collisionCheck)).second;
  m_fwdSamples.push_back(m_fwdRoot);
  m_fwdFront.push_back(m_fwdRoot);

  //closest configs and distance measures
  m_minDistance = dist_to_objective(m_fwdRoot);
  m_closestFwdSample = m_fwdRoot;
  m_addedToFront = true;
  m_numSamples = 0;
}

DEERPlanner::~DEERPlanner() {
  Configuration *pSmp;
  for (auto &m_fwdSample : m_fwdSamples) {
    pSmp = m_fwdSample;
    delete pSmp;
  }
}

std::list<Configuration *> &DEERPlanner::getSamples() {
  return m_fwdSamples;
}


void DEERPlanner::generateSamples() {

  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  bool readyToFix = false;

  /// %%%%%%%%%%%%%%%%%%%%% First downhill approach as far as possible
  log("samplingStatus") << "Starting downhill moves"<<endl;
  approachTarget();

  double end_time = timer.ElapsedTime();

  if ( m_minDistance < m_convergeDistance){
    log("samplingStatus") << "Reached DEER distances down to "<<m_minDistance<<" in "<<end_time - start_time<<" seconds. "<<endl;
    readyToFix = true;
  }
  /// %%%%%%%%%%%%%%%%%%%%%

  start_time = timer.LastElapsedTime();

  /// %%%%%%%%%%%%%%%%%%%%% Use motion planning with moving front to approach further
  log("samplingStatus")<<endl;
  log("samplingStatus") << "Starting motion planning"<<endl;
  navigateToTarget();

  end_time = timer.ElapsedTime();

  if(m_minDistance <= m_convergeDistance) {
    //Set additional distance constraint
    log("samplingStatus") << "Reached DEER distances down to " << m_minDistance << " in " << end_time - start_time
                          << " seconds. " << endl;
    readyToFix = true;
  }

  if( readyToFix ){
    /// %%%%%%%%%%%%%%%%%%%%% Use randomized exploration with new DEER distance constraints
//    log("samplingStatus") << "To implement: Fixing DEER constraints and exploring"<<endl;

//    double startTime = timer.LastElapsedTime();
//
//    fixConstraints();
//    exploreAroundTarget();
//
//    double endTime = timer.ElapsedTime();
//
//    log("samplingStatus") << "Finished exploring in "<<endTime - startTime<<" seconds. "<<endl;
    /// %%%%%%%%%%%%%%%%%%%%%
  }
  else{
    /// %%%%%%%%%%%%%%%%%%%%% Couldn't reach target distances, return
    log("samplingStatus") << "Did not reach Deer distances to below "<<m_minDistance<<" within "<<m_stopAfter
    <<" samples. Closest conformation "<<m_closestFwdSample->m_id<<" at distance "<<m_minDistance<<endl;
  }
}

void DEERPlanner::approachTarget() {

  gsl_vector* gradient = gsl_vector_alloc(m_protein->m_spanningTree->getNumDOFs());

  double previousDist = dist_to_objective(m_fwdSamples.back());
  double dist;
  while (m_numSamples < m_stopAfter) {

    if (previousDist<m_convergeDistance) break;

    Configuration* seed = m_fwdSamples.back();
    m_directedDirection->gradient(seed, nullptr, gradient); //directed move

    Configuration* newConf = m_directedMove->move(seed, gradient);
    if (newConf->updatedMolecule()->inCollision()) {
      cout<<"Could not find new conformation"<<endl;
      delete newConf;
      break; /// end of downhill approach
    }
    else {
      dist = dist_to_objective(newConf);
      if(dist > 0.999*previousDist){
        log("samplingStatus")<<"Not getting closer, ending downhill approach"<<endl;
        delete newConf;
        break; /// end of downhill approach
      }
      else {
        m_numSamples++;
        newConf->m_id = m_numSamples;
        m_fwdSamples.push_back(newConf);

        double violationNorm = m_protein->checkCycleClosure(newConf);

        evaluateDistances(newConf);
        updateFwdFront(newConf);

        log("samplingStatus") << "> New structure: " << newConf->getMolecule()->getName() << "_new_" << m_numSamples << ".pdb";
        log("samplingStatus") << " .. Dist initial: " << setprecision(6) << newConf->m_distanceToIni;
        log("samplingStatus") << " .. Distance to objective: " << setprecision(6) << newConf->m_distanceToTarget;
        log("samplingStatus") << " .. access dofs: " << newConf->m_clashFreeDofs;
        log("samplingStatus") << " .. norm constr. viol.: " << violationNorm << endl;
        IO::writeNewSample(newConf, m_fwdRoot, newConf->m_id, m_workingDir, m_saveData);

        previousDist = dist;
      }
    }
  }
}

void DEERPlanner::navigateToTarget() {

  static int failedTrials = 0, totalTrials = 0;

  Configuration *qTarget = nullptr, *qSeed = nullptr, *qNew = nullptr; //this qTarget is either global or random and can change for each sample

  while (m_numSamples < m_stopAfter) {//sample at most until the number of samples has been reached
    ++totalTrials;

    CTKTimer timer;
    double start_time = timer.getTimeNow();

    double randVal = Random01(); // used to determine bias to target

//    if (qNew == nullptr || m_addedToFront == false || randVal > m_biasToTarget ) {
    if (qNew == nullptr || randVal > m_biasToTarget){//create new target and compute closest seed

      if (qTarget != nullptr) {//prevent leakage, delete existing target configuration from previous round
        delete qTarget;
      }

      qTarget = GenerateRandConf();//New random target configuration
      qSeed = SelectSeed(qTarget); //Compute closest seed
    }
    else{//directed move towards target
      qSeed = m_closestFwdSample; //qSeed = m_fwdSamples.back();
    }
    double end_time = timer.getTimeNow();
    selectNodeTime += end_time - start_time;


    gsl_vector *gradient = gsl_vector_calloc(m_protein->m_spanningTree->getNumDOFs());
//    log("planner") << "Base dofs " << qSeed->m_dofs <<endl;
    if (m_isBlended) {
      auto &blendedDir = reinterpret_cast<BlendedDirection &>(*m_direction);
      blendedDir.changeWeight(0, double(m_numSamples) / double(m_stopAfter));
      blendedDir.changeWeight(1, 1.0 - double(m_numSamples) / double(m_stopAfter));
    }

    ///Two types of move
    if(randVal < m_biasToTarget) {//directed move
      m_directedDirection->gradient(qSeed, nullptr,gradient);
      qNew = m_directedMove->move(qSeed, gradient); //Perform move
    }
    else{//random move
      m_direction->gradient(qSeed, qTarget, gradient); //computes the search m_direction for a new sample
//      gsl_vector_scale_to_length(gradient, m_maxRotation);
      qNew = m_move->move(qSeed, gradient); //Perform move
    }

    if (qNew->updatedMolecule()->inCollision()) {
      failedTrials++;
      delete qNew;
      qNew = nullptr;

    } else {//collision-free
      m_numSamples++;
      gsl_vector_free(gradient);

      //Push-back in forward tree
      m_fwdSamples.push_back(qNew);
      qNew->m_id = m_numSamples;

      //Potentially reject new config if large violations?
      double violationNorm = m_protein->checkCycleClosure(qNew);

      //Distance computations
      evaluateDistances(qNew);

      //Enthalpy and entropy computation
//      log("planner")<<"Length of all collisions: "<<allCollisions.size()<<endl;
//      pair<double,double> enthalpyVals=m_molecule->vdwEnergy(&allCollisions,m_options.collisionCheck);
//      qNew->m_vdw = enthalpyVals.first;
//      qNew->m_deltaH = enthalpyVals.second - new_q->m_vdw;

      qNew->m_vdwEnergy = qNew->getMolecule()->vdwEnergy(m_collisionCheck);

      log("planner") << "> New structure: " << qNew->getMolecule()->getName() << "_new_" << m_numSamples
                     << ".pdb, accessible dofs: " << qNew->m_clashFreeDofs << endl << endl;

      log("samplingStatus") << "> New structure: " << qNew->getMolecule()->getName() << "_new_" << m_numSamples << ".pdb";
      log("samplingStatus") << " .. Dist initial: " << setprecision(6) << qNew->m_distanceToIni;
      log("samplingStatus") << " .. Dist target: " << setprecision(6) << qNew->m_paretoFrontDistance;
      log("samplingStatus") << " .. Closest distance: " << setprecision(6) << m_minDistance;
      log("samplingStatus") << " .. access dofs: " << qNew->m_clashFreeDofs;
      log("samplingStatus") << " .. norm constr. viol.: " << violationNorm << endl;

      IO::writeNewSample(qNew, m_fwdRoot, qNew->m_id, m_workingDir, m_saveData);

      //Check if front needs to be updated, deletes nullspace if not necessary anymore
      updateFwdFront(qNew);

      //Check if target has been reached!
      if (m_minDistance < m_convergeDistance) {
        log() << "Reached target, not creating more samples!" << " distance: " << m_minDistance << endl;
        break;
      }
    }
  }

  log() << endl << "Max tree depth (excluding the m_root) = " << m_max_depth << endl;
  log() << "failed " << failedTrials << " times due to collision." << endl;

}

void DEERPlanner::fixConstraints() {
//  /// Generate new molecule with DEER distance constraints
//  m_closestFwdSample->updatedMolecule();
//  Molecule* exploreProtein = m_protein->deepClone();
//
//  double d=0;
//  for (auto &m_goalDistance : m_goalDistances) {
//    Atom *a1 = exploreProtein->getAtom(get<0>(m_goalDistance)->getId());
//    Atom *a2 = exploreProtein->getAtom(get<1>(m_goalDistance)->getId());
//    DBond *new_db = new DBond(a1, a2);
//    exploreProtein->addDBond(new_db);
//  }
//
//  Selection resNetwork("all");
//  exploreProtein->initializeTree(resNetwork);
}

void DEERPlanner::exploreAroundTarget() {

}

void DEERPlanner::updateFwdFront(Configuration *qNew) {

  std::list<Configuration *>::iterator cit;
  int i = 1;
  m_addedToFront = false;

  for (cit = m_fwdFront.begin(); cit != m_fwdFront.end(); cit++) {
    if (qNew->m_paretoFrontDistance <= (*cit)->m_paretoFrontDistance) {
//      log("planner") << "Inserting into pareto front at: " << i << ", with distance: " << qNew->m_paretoFrontDistance << endl;
      m_fwdFront.insert(cit, qNew); //potentially use it as new seed / target for reverse planner
      m_addedToFront = true;
      break;
    }
    i++;
  }
  if (cit == m_fwdFront.end() && m_fwdFront.size() < m_frontSize) {
    m_fwdFront.push_back(qNew);
//    log("planner") << "Inserting at the end at: " << i << ", with distance: " << qNew->m_paretoFrontDistance << endl;
    m_addedToFront = true;
  }
  if (m_fwdFront.size() > m_frontSize)
    m_fwdFront.pop_back();//keep length at maximum

  ///Configurations outside of front will never be used as seeds, delete nullspace
  if(!m_addedToFront){
    qNew->deleteNullspace();
  }
}

void DEERPlanner::evaluateDistances(Configuration *qNew) {

  ///use this if real target structures are supplied
  double dist = dist_to_objective(qNew);

  qNew->m_distanceToTarget = dist;

  //Distance to own m_root
  qNew->m_distanceToIni = m_metric->distance(qNew, m_fwdRoot);

  //Distance to parents
  if (qNew->getParent() != nullptr) {
    qNew->m_distanceToParent = m_metric->distance(qNew, qNew->getParent());
  }
///adapt this to closest reverse structure if supplied
  qNew->m_paretoFrontDistance = dist;
  if (qNew->m_paretoFrontDistance < m_minDistance) {//update current shortest configs
    m_minDistance = qNew->m_paretoFrontDistance;
    m_closestFwdSample = qNew;
  }

  if (qNew->m_treeDepth > m_max_depth)
    m_max_depth = qNew->m_treeDepth;
}

Configuration * DEERPlanner::GenerateRandConf() {

  auto *pNewSmp = new Configuration(m_protein);
  size_t numDofs = m_protein->m_spanningTree->getNumDOFs();
//  double length = 0.0;
  for (int i = 0; i < numDofs; ++i) {
    pNewSmp->m_dofs[i] = Math3D::dPi * RandomN1P1();
//    length += pNewSmp->m_dofs[i] * pNewSmp->m_dofs[i];
  }
//  length = sqrt(length);

  pNewSmp->m_id = -1;
  return pNewSmp;
}

Configuration *DEERPlanner::SelectSeed(Configuration *pTarget) {

  Configuration *pMinSmp;
  double minDistance = 1000000.0;
  double distance;
  pMinSmp = m_fwdFront.back(); //prevent hole in case distance greater than minDistance

  pTarget->updatedMolecule();//update target protein

  for (auto const &pSmp : m_fwdFront) {
    distance = m_metric->distance(pSmp, pTarget);
    if (distance < minDistance) {
      minDistance = distance;
      pMinSmp = pSmp;
    }
  }
  return pMinSmp;
}

void DEERPlanner::createTrajectory() {

  m_closestFwdSample->updateMolecule();

  log() << "The forward trajectory ends at sample " << m_closestFwdSample->m_id << endl;

  const string &out_path = m_workingDir;
  const string &name = m_protein->getName();

  string out_collPdb = out_path + "output/" + name + "_path.pdb";
  ///save pyMol movie script
  string out_pyMol = out_path + "output/" + name + "_pathMov.pml";

  IO::writeTrajectory(m_protein, out_collPdb, out_pyMol, nullptr);
}

double DEERPlanner::dist_to_objective(Configuration* conf)
{
  conf->updatedMolecule();
  double d=0;
  for (auto &m_goalDistance : m_goalDistances) {
    Coordinate c1= get<0>(m_goalDistance)->m_position;
    Coordinate c2= get<1>(m_goalDistance)->m_position;

//    d=d+fabs((sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z))-get<2>(m_goalDistances[i])));
    double dterm = sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z))-get<2>(m_goalDistance);
    d += dterm*dterm;
  }
//  return d;
  return sqrt(d);
}