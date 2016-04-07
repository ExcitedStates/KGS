//
// Created by Dominik Budday on 09.03.16.
//

#ifndef KGS_BIDIRECTIONALMOVINGFRONT_H
#define KGS_BIDIRECTIONALMOVINGFRONT_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "metrics/Metric.h"
#include "SamplingOptions.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"
#include "directions/Direction.h"

class Molecule;
class Configuration;


/**
 * A sampling planner based on a bidirectional RRT (paper KuffnerLaValle 2000), adapted with an iteratively
 * updated moving front from which samples and seeds are selected efficiently. The planner is normally directed,
 * sampling from a provided initial and target structure. We generate a random target / select a target from the
 * reverse front, use the closest sample in the forward front as seed and move towards the target. If our new sample
 * is a) in the moving front, we select it as new seed and continue the search.
 *    b) not in the moving front, we generate a new target and seed to continue in the same direction
 *    c) not accepted (clashing), we switch search directions.
 * The planner continues until a given distance threshold between the fronts is reached, or a given number of samples.
 */

class BidirectionalMovingFront : public SamplingPlanner{
 public:
  BidirectionalMovingFront(Molecule * protein, Move& move, Direction& direction, Molecule * target);
  ~BidirectionalMovingFront();

  void GenerateSamples();

  ConfigurationList& Samples(){ return m_fwdSamples; }

  ConfigurationList& getFwdSamples(){ return m_fwdSamples; }
  ConfigurationList& getRevSamples(){ return m_revSamples; }
  ConfigurationList& getFwdFront(){ return m_fwdFront; }
  ConfigurationList& getRevFront(){ return m_revFront; }

  void createTrajectory(); ///< overwrites the parent function createTrajectory

 protected:
  Configuration* GenerateRandConf();
  Configuration* SelectSeed(Configuration *pTarget);

  void setClosestConfigRMSD();

  //metrics::Metric* m_metric;
  metrics::Metric* rmsdMetric;

  Direction& direction;

 private:

  const int stopAfter;

  void updateFwdFront(Configuration* qNew);
  void evaluateDistances(Configuration* qNew);
  void swapFwdRev();

  Molecule *m_protein;
  Molecule *m_target;

  metrics::Metric* m_metric;            ///< Metric used to measure distances between samples

  Configuration* m_fwdRoot;
  Configuration* m_revRoot;

  ConfigurationList m_fwdSamples;     ///< all forward samples
  ConfigurationList m_revSamples;   ///< all reverse samples
  ConfigurationList m_fwdFront;      ///< forward moving front
  ConfigurationList m_revFront;     ///< reverse moving front

  Configuration* m_currentGlobalTarget; // this is the target for this round
  Configuration* m_closestFwdSample; // own configuration
  Configuration* m_closestRevSample; //closest target configuration

  int m_sample_id;
  int m_max_depth;
  double m_minDistance;

  bool m_addedToFront;
  int m_frontSize;

  int m_nCDCall = 0;
  int m_nMetricsCall = 0;

  int movesRejected = 0;
  int movesAccepted = 0;
};


#endif //KGS_BIDIRECTIONALMOVINGFRONT_H
