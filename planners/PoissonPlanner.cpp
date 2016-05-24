/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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

#include "PoissonPlanner.h"

#include <iomanip>
#include <stack>
#include <directions/Direction.h>
#include <directions/RandomDirection.h>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "Logger.h"
#include "core/Transformation.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "metrics/Dihedral.h"

using namespace std;

PoissonPlanner::PoissonPlanner(Molecule * protein, Move& move, metrics::Metric& metric):
	SamplingPlanner(move,metric),
  stop_after(SamplingOptions::getOptions()->samplesToGenerate),
  max_rejects_before_close(SamplingOptions::getOptions()->poisson_max_rejects_before_close),
  m_bigRad(SamplingOptions::getOptions()->stepSize*4.0/3.0),
  m_lilRad(m_bigRad/2),
  protein(protein)
{
  m_root = new Configuration( protein );
  m_root->updateMolecule();
  m_root->computeCycleJacobianAndNullSpace();
  m_root->m_id = 0;
  open_samples.push_back( m_root );
  all_samples.push_back( m_root );
  updateMaxDists(m_root,m_root);
}

PoissonPlanner::~PoissonPlanner() {
  Configuration *pSmp;
  for (auto iter=all_samples.begin(); iter!=all_samples.end(); iter++) {
    pSmp = *iter;
    delete pSmp;
  }
}



void PoissonPlanner::GenerateSamples()
{
  //cout<<"PoissonPlanner::GenerateSamples()"<<endl;
  Direction* direction = new RandomDirection();
  gsl_vector* gradient = gsl_vector_alloc(protein->totalDofNum());
  double origStepSize = move.getStepSize();

  int sample_num = 0;
  int rejected_clash     = 0;
  int rejected_collision = 0;

  while(sample_num<stop_after && !open_samples.empty()) {
    //Pick random open conformation
    auto it = open_samples.begin();
    std::advance(it, rand()%open_samples.size());
    Configuration* seed = *it;

    vector<Configuration*> nearSeed;
    if(!m_checkAll) {
      collectPossibleChildCollisions(seed, nearSeed, m_root);
    }


    //Make max_rejects_before_close attempts at perturbing it
    size_t attempt;
    for( attempt=0; attempt<max_rejects_before_close; attempt++ ) {
      //cout<<"PoissonPlanner::GenerateSamples() - attempt "<<attempt<<endl;
      move.setStepSize(origStepSize);
      direction->gradient(seed, nullptr, gradient); // Compute random gradient
      Configuration *pert = move.move(seed, gradient); //Perform move

      // Scale gradient so move is in Poisson disc
      double dist = m_metric.distance(pert, seed);
      //log("samplingStatus")<<" - seed-to-new distance is now "<<dist<<" (should be between "<<m_lilRad<<" and "<<m_bigRad<<")"<<endl;
      //cout<<"PoissonPlanner - dist: "<<dist<<", lil: "<<m_lilRad<<", big: "<<m_bigRad<<endl;
      int scaleAttempts = 0;
      while( dist<m_lilRad || dist>m_bigRad){
        if(++scaleAttempts==5) break;
        double gradientScale = (m_bigRad+m_lilRad)/(2.0*dist);
        //cout<<"PoissonPlanner - scaling up/down by "<<gradientScale<<endl;
        gsl_vector_scale(gradient, gradientScale);
        //move.setStepSize(move.getStepSize()*gradientScale);
        delete pert;
        pert = move.move(seed, gradient);
        dist = m_metric.distance(pert, seed);
        //cout<<"PoissonPlanner -  dist is now: "<<dist<<endl;
        //log("samplingStatus")<<" - seed-to-new distance is now "<<dist<<" (should be between "<<m_lilRad<<" and "<<m_bigRad<<")"<<endl;
      }

      if(scaleAttempts==5){
        delete pert;
        continue;
      }

      //If clashing just continue
      if(pert->updatedMolecule()->inCollision() ) {
        rejected_clash++;
        delete pert;
        //log("samplingStatus")<<" - rejected from clash)"<<endl;
        continue;
      }


      //Check if close to existing
      bool too_close_to_existing = false;
      if(m_checkAll) {
        for (auto const &v: all_samples) {
          if (m_metric.distance(pert, v) < m_lilRad) {
            too_close_to_existing = true;
            break;
          }
        }
      }else{
        for (auto const &v: nearSeed) {
          if (m_metric.distance(pert, v) < m_lilRad) {
            too_close_to_existing = true;
            break;
          }
        }
      }
      if (too_close_to_existing) {
        rejected_collision++;
        delete pert;
        //log("samplingStatus")<<" - rejected from collision"<<endl;
        continue;
      }

      //Perturbation is OK: Add to open
      pert->m_id = sample_num;
      open_samples.push_back(pert);
      all_samples.push_back(pert);
      writeNewSample(pert, all_samples.front(), sample_num);
      pert->m_distanceToIni    = m_metric.distance(pert,all_samples.front());
      pert->m_distanceToParent = m_metric.distance(pert,seed);

      log("samplingStatus") << "> New structure: "<<pert->getMolecule()->getName()<<"_new_"<<sample_num<<".pdb";
      log("samplingStatus") << " .. initial dist.: "<< setprecision(6)<<pert->m_distanceToIni;
      log("samplingStatus") << endl;

      sample_num++;
      break;
    }

    //If loop is not terminated with a break it means the seed should be closed
    if(attempt==max_rejects_before_close){
      open_samples.erase(it);
      closed_samples.push_back(seed);
    }
  }
  log("samplingStatus")<<"Poisson-planner: "<<open_samples.size()<<" open samples and ";
  log("samplingStatus")<<closed_samples.size()<<" closed samples on termination"<<endl;
  log("samplingStatus")<<"Poisson-planner: Rejects from clash: "<<rejected_clash<<endl;
  log("samplingStatus")<<"Poisson-planner: Rejects from tree-collision: "<<rejected_collision<<endl;

  gsl_vector_free(gradient);
  delete direction;

}

double PoissonPlanner::memo_distance(Configuration* c1, Configuration* c2)
{
  if(c1==c2) return 0.0;

  auto it = m_distances.find({c1,c2});
  if( it!=m_distances.end() ) return it->second;

  double dist = m_metric.distance(c1,c2); //Expensive
  m_distances[{c1,c2}] = dist;
  m_distances[{c2,c1}] = dist;

  return dist;
}

void PoissonPlanner::collectPossibleChildCollisions(
    Configuration* conf,
    std::vector<Configuration*>& ret
)
{
  collectPossibleChildCollisions(conf, ret, m_root);
}

void PoissonPlanner::collectPossibleChildCollisions(
    Configuration* conf,
    std::vector<Configuration*>& ret,
    Configuration* v
)
{
  double d = memo_distance(v,conf);
  if( d >= m_maxDist[v]+m_bigRad+m_lilRad )
    return;

  if(d<m_bigRad+m_lilRad) ret.push_back(v);

  for(auto const& child: v->getChildren())
    collectPossibleChildCollisions(conf, ret, child);

}

void PoissonPlanner::updateMaxDists(Configuration* newConf)
{
  updateMaxDists(newConf, newConf);
}
void PoissonPlanner::updateMaxDists(Configuration* v, Configuration* newConf)
{
  double d = memo_distance(v,newConf);
  if( d>m_maxDist[v] ) m_maxDist[v] = d;

  if(v->getParent())
    updateMaxDists(v->getParent(), newConf);
}
