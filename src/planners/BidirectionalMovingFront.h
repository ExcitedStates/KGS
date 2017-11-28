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
// Created by Dominik Budday on 09.03.16.
//

#ifndef KGS_BIDIRECTIONALMOVINGFRONT_H
#define KGS_BIDIRECTIONALMOVINGFRONT_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <metrics/RMSD.h>

#include "metrics/Metric.h"
#include "core/Molecule.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"
#include "directions/Direction.h"


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
  BidirectionalMovingFront(
      Molecule *protein,
      Direction *direction,
      Molecule *target,
      Selection &metricSelection,
      std::string& collisionCheck,
      bool blendedDir,
      int stopAfter,
      int frontSize,
      double stepSize,
      int switchAfter,
      double convergeDistance,
      bool alignAlways,
      double biasToTarget
  );

  ~BidirectionalMovingFront();

  void generateSamples();

  std::list<Configuration*>& getSamples();

  std::list<Configuration*>& getFwdSamples(){ return m_fwdSamples; }
  std::list<Configuration*>& getRevSamples(){ return m_revSamples; }
  std::list<Configuration*>& getFwdFront(){ return m_fwdFront; }
  std::list<Configuration*>& getRevFront(){ return m_revFront; }

  void createTrajectory(); ///< overwrites the parent function createTrajectory

 protected:
  Configuration* GenerateRandConf();
  Configuration* SelectSeed(Configuration *pTarget);

  void setClosestConfigRMSD();

  Direction* direction;

 private:

  const int m_stopAfter;

  Molecule *m_protein;
  Molecule *m_target;

  Configuration* m_fwdRoot;
  Configuration* m_revRoot;

  std::list<Configuration*> m_fwdSamples;     ///< all forward samples
  std::list<Configuration*> m_revSamples;   ///< all reverse samples
  std::list<Configuration*> m_fwdFront;      ///< forward moving front
  std::list<Configuration*> m_revFront;     ///< reverse moving front

  Configuration* m_currentGlobalTarget; // this is the target for this round
  Configuration* m_closestFwdSample; // own configuration
  Configuration* m_closestRevSample; //closest target configuration

  metrics::RMSD* m_rmsd; ///< Used for aligning

  int m_max_depth;
  double m_minDistance;

  bool m_addedToFront;
  int m_frontSize;
  bool m_isBlended;

  std::string m_collisionCheck;
  bool m_samplingForward;
  double m_stepSize;
  int m_switchAfter;
  double m_convergeDistance;
  bool m_alignAlways;
  double m_biasToTarget;

  int m_nCDCall = 0;
  int m_nMetricsCall = 0;

  int movesRejected = 0;
  int movesAccepted = 0;

  void updateFwdFront(Configuration* qNew);
  void evaluateDistances(Configuration* qNew);
  void swapFwdRev();

};


#endif //KGS_BIDIRECTIONALMOVINGFRONT_H
