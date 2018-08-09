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

#ifndef KGS_DEERPLANNER_H
#define KGS_DEERPLANNER_H


#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "metrics/Metric.h"
#include "core/Molecule.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"
#include "directions/Direction.h"


/**
 * A sampling planner based on the RRT, adapted with an iteratively
 * updated moving front from which samples and seeds are selected efficiently. The planner is normally directed,
 * sampling from a provided initial structure. The DEER distance distributions are used as input to compute
 * a set of target distances or randomized distances, use the closest sample in the forward front as seed and move towards the target. 
 * If our new sample is
 *    a) in the moving front, we select it as new seed and continue the search.
 *    b) not in the moving front, we generate a new target and seed to continue in the same direction
 *    c) not accepted (clashing), we do something.
 * The planner continues until a given distance threshold between the fronts is reached, or a given number of samples.
 */


class DEERPlanner : public SamplingPlanner{
 public:
  DEERPlanner(
      Molecule *protein,
      Direction *direction,
      Move *directedMove,
      Direction *directedDirection,
      std::vector< std::tuple<Atom*, Atom*, double> > goal_distances,
      std::string& collisionCheck,
      bool blendedDir,
      int stopAfter,
      int frontSize,
      double stepSize,
      double convergeDistance,
      double biasToTarget
  );

  ~DEERPlanner();

  void generateSamples();

  std::list<Configuration*>& getSamples();

  std::list<Configuration*>& getFwdSamples(){ return m_fwdSamples; }
  std::list<Configuration*>& getFwdFront(){ return m_fwdFront; }

  Configuration* getClosestConfig(){return m_closestFwdSample; }

  void createTrajectory(); ///< overwrites the parent function createTrajectory

 protected:
  Configuration* GenerateRandConf();
  Configuration* SelectSeed(Configuration *pTarget);

  Direction* m_direction;

  void approachTarget();
  void navigateToTarget();
  void exploreAroundTarget();
  void fixConstraints();

 private:

  const int m_stopAfter;

  Molecule *m_protein;

  Configuration* m_fwdRoot;

  std::vector< std::tuple<Atom*, Atom*, double> > m_goalDistances;

  std::list<Configuration*> m_fwdSamples;     ///< all forward samples
  std::list<Configuration*> m_fwdFront;      ///< forward moving front

  Configuration* m_closestFwdSample; // own configuration

  Move *m_directedMove;
  Direction *m_directedDirection;

  int m_max_depth;
  double m_minDistance;

  bool m_addedToFront;
  int m_frontSize;
  bool m_isBlended;
  int m_numSamples;

  std::string m_collisionCheck;
  double m_stepSize;
  double m_convergeDistance;
  double m_biasToTarget;

  void updateFwdFront(Configuration* qNew);
  void evaluateDistances(Configuration* qNew);

  double dist_to_objective(Configuration* conf);

};



#endif //KGS_DEERPLANNER_H
