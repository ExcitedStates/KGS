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



#include "PoissonPlanner2.h"

#include <iomanip>
#include <stack>
#include <directions/Direction.h>
#include <directions/RandomDirection.h>
#include <loopclosure/ExactIK.h>
#include <IO.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "Logger.h"
#include "core/Transformation.h"

using namespace std;

PoissonPlanner2::PoissonPlanner2(
    Molecule * protein,
    vector<ResTriple>& exactIKtriples,
    int stopAfter,
    int maxRejects,
    double stepSize,
    const string &resNetwork
):
    SamplingPlanner(),
//    m_stopAfter(ExploreOptions::getOptions()->samplesToGenerate),
//    m_maxRejectsBeforeClose(ExploreOptions::getOptions()->poissonMaxRejectsBeforeClose),
//    m_bigRad(ExploreOptions::getOptions()->stepSize*4.0/3.0),
    m_stopAfter(stopAfter),
    m_maxRejectsBeforeClose(maxRejects),
    m_bigRad(stepSize*4.0/3.0),
    m_lilRad(m_bigRad/2),
    m_protein(protein),
    m_resNetwork(resNetwork),
    m_ikTriples(exactIKtriples)
{
  m_root = new Configuration( protein );
  m_root->m_id = 0;
  open_samples.push_back( m_root );
  all_samples.push_back( m_root );
  updateMaxDists(m_root);
}

PoissonPlanner2::~PoissonPlanner2() {
  Configuration *pSmp;
  for (auto iter=all_samples.begin(); iter!=all_samples.end(); iter++) {
    pSmp = *iter;
    delete pSmp;
  }
}



void PoissonPlanner2::generateSamples()
{
  //cout<<"PoissonPlanner2::generateSamples()"<<endl;
  Selection sel(m_resNetwork);
  Direction* direction = new RandomDirection(sel);
  gsl_vector* gradient = gsl_vector_alloc(m_protein->totalDofNum());
//  double oriMaxRotation = m_move->getMaxRotation();
  m_move->setScalingFlag(false);

  int sample_num = 1;
  int rejected_clash     = 0;
  int rejected_collision = 0;

  while(sample_num<m_stopAfter && !open_samples.empty()) {
    //Pick random open conformation
    auto it = open_samples.begin();
    std::advance(it, rand()%open_samples.size());
    Configuration* seed = *it;
    log("samplingStatus") << "Using configuration "<<seed->m_id<<" as seed. "<<open_samples.size()<<" open, "<<closed_samples.size()<<" closed samples"<<endl;

    vector<Configuration*> nearSeed;
    if(!m_checkAll) {
      collectPossibleChildCollisions(seed, nearSeed, m_root, m_bigRad);
    }


    //Make m_maxRejectsBeforeClose attempts at perturbing it
    size_t attempt;
    for( attempt=0; attempt<m_maxRejectsBeforeClose; attempt++ ) {
//      m_move->setMaxRotation(oriMaxRotation);
      direction->gradient(seed, nullptr, gradient); // Compute random gradient
      Configuration *pert = m_move->move(seed, gradient); //Perform move

      // Scale gradient so move is in Poisson disc
      double dist = m_metric->distance(pert, seed);
      int scaleAttempts = 0;
      while( dist<m_lilRad || dist>m_bigRad){
        if(++scaleAttempts==5) break;
        double gradientScale = (m_bigRad+m_lilRad)/(2.0*dist);
        gsl_vector_scale(gradient, gradientScale);
        delete pert;
        pert = m_move->move(seed, gradient);
        dist = m_metric->distance(pert, seed);
      }

      if(scaleAttempts==5){
        delete pert;
        continue;
      }

      //If clashing just continue
      if(pert->updatedMolecule()->inCollision() ) {
        rejected_clash++;
        delete pert;
        continue;
      }


      //Check if close to existing
      bool too_close_to_existing = false;
      for (auto const &v: nearSeed) {
        double dist = memo_distance(pert, v);
        if (dist < m_lilRad) {
          too_close_to_existing = true;
          break;
        }
      }
      if (too_close_to_existing) {
        rejected_collision++;
        delete pert;
        continue;
      }

      //Perturbation is OK: Add to open
      pert->m_id = sample_num;
      open_samples.push_back(pert);
      all_samples.push_back(pert);
      pert->m_distanceToIni    = memo_distance(pert,m_root);
      pert->m_distanceToParent = memo_distance(pert,seed);
      updateMaxDists(pert);
      IO::writeNewSample(pert, m_root, sample_num, m_workingDir, m_saveData);

      log("samplingStatus") << "> "<<pert->getMolecule()->getName()<<"_new_"<<sample_num<<".pdb";
      log("samplingStatus") << " .. init dist: "<< setprecision(3)<<pert->m_distanceToIni;
      log("samplingStatus") << " .. seed dist: "<< setprecision(3)<<pert->m_distanceToParent;
      log("samplingStatus") << endl;

      sample_num++;
    }

    //Perform exact IK rebuild
    ExactIK ik;
    for(auto triple: m_ikTriples) {
      seed->updateMolecule();
      vector<Configuration *> rebuilt = ik.rebuildLoop(get<0>(triple), get<1>(triple), get<2>(triple) );
      for(Configuration* pert: rebuilt){

        //Collision with seed is quick to check
        pert->m_distanceToParent = memo_distance(pert,seed);
        if(pert->m_distanceToParent<m_lilRad){
          rejected_collision++;
          delete pert;
          continue;
        }

        //If clashing just continue
        if(pert->updatedMolecule()->inCollision() ) {
          rejected_clash++;
          delete pert;
          continue;
        }

        //Check if close to existing
        vector<Configuration*> nearPert;
        collectPossibleChildCollisions(pert, nearPert, 0.0);
        bool too_close_to_existing = false;
        for (auto const &v: nearPert) {
//        for (auto const &v: all_samples) {
          double dist = memo_distance(pert, v);
          if (dist < m_lilRad) {
            too_close_to_existing = true;
            break;
          }
        }
        if (too_close_to_existing) {
          rejected_collision++;
          delete pert;
          continue;
        }

        //Perturbation is OK: Add to open
        pert->m_id = sample_num;
        pert->m_distanceToIni    = memo_distance(pert,m_root);
        open_samples.push_back(pert);
        all_samples.push_back(pert);
        updateMaxDists(pert);
        IO::writeNewSample(pert, m_root, sample_num, m_workingDir, m_saveData);

        log("samplingStatus") << "> "<<pert->getMolecule()->getName()<<"_new_"<<sample_num<<".pdb";
        log("samplingStatus") << " .. init dist: "<< setprecision(3)<<pert->m_distanceToIni;
        log("samplingStatus") << " .. seed dist: "<< setprecision(3)<<pert->m_distanceToParent;
        log("samplingStatus") << " .. from exact IK";
        log("samplingStatus") << endl;

        sample_num++;

      }

    }

    open_samples.erase(it);
    closed_samples.push_back(seed);
  }
  log("samplingStatus")<<"Poisson-planner: "<<open_samples.size()<<" open samples and ";
  log("samplingStatus")<<closed_samples.size()<<" closed samples on termination ";
  log("samplingStatus")<<"("<<all_samples.size()<<" total)"<<endl;
  log("samplingStatus")<<"Poisson-planner: Rejects from clash:          "<<rejected_clash<<endl;
  log("samplingStatus")<<"Poisson-planner: Rejects from tree-collision: "<<rejected_collision<<endl;

  gsl_vector_free(gradient);
  delete direction;

}

double PoissonPlanner2::memo_distance(Configuration* c1, Configuration* c2)
{
  if(c1==c2) return 0.0;

  auto it = m_distances.find({c1,c2});
  if( it!=m_distances.end() ) return it->second;

  double dist = m_metric->distance(c1,c2); //Expensive
  m_distances[{c1,c2}] = dist;
  m_distances[{c2,c1}] = dist;

  return dist;
}

void PoissonPlanner2::collectPossibleChildCollisions(
    Configuration* conf,
    std::vector<Configuration*>& ret,
    double childOffset
)
{
  collectPossibleChildCollisions(conf, ret, m_root, childOffset);
}

void PoissonPlanner2::collectPossibleChildCollisions(
    Configuration* conf,
    std::vector<Configuration*>& ret,
    Configuration* v,
    double childOffset
)
{
  if(v->m_id<0) return;

  double d = memo_distance(v,conf);
  //if( d >= m_maxDist[v]+m_bigRad+m_lilRad )
  if( d >= m_maxDist[v]+childOffset+m_lilRad )
    return;

  //if(d<m_bigRad+m_lilRad) ret.push_back(v);
  if(d<childOffset+m_lilRad) ret.push_back(v);

  for(auto const& child: v->getChildren()) {
//    cout<<"collectPossibleChildCollisions - "<<child->m_id<<" is a child of "<<v->m_id<<endl;
    collectPossibleChildCollisions(conf, ret, child, childOffset);
  }

}

void PoissonPlanner2::updateMaxDists(Configuration* newConf)
{
  updateMaxDists(newConf, newConf);
}
void PoissonPlanner2::updateMaxDists(Configuration* v, Configuration* newConf)
{
  double d = memo_distance(v,newConf);
  if( d>m_maxDist[v] ) m_maxDist[v] = d;

  Configuration* parent = v->getParent();
  if(parent!=nullptr) updateMaxDists(parent, newConf);
}
