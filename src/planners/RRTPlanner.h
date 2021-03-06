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

#ifndef RRTPLANNER_H
#define RRTPLANNER_H

#include <iostream>
#include <fstream>
#include <string>
#include <directions/Direction.h>

#include "metrics/Metric.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"


class Molecule;

#define MAX_BUCKET_NUM 101
#define NUM_BINS 100
#define DEFAULT_MAX_RMSD 25


class RRTPlanner : public SamplingPlanner {
 public:
  RRTPlanner(
      Molecule *mol,
      Direction *direction,
      double radius,
      int numSamples,
      int gradientSelection,
      double stepSize,
      double maxRotation,
      bool scaleToRadius
  );

  ~RRTPlanner();

  void generateSamples();

  std::list<Configuration *> &getSamples() { return m_samples; }

 protected:
  Configuration *GenerateRandConf();

  Configuration *SelectNodeFromBuckets(Configuration *pTarget);

  unsigned int m_current_max_bucket_id;
  std::list<Configuration *> distance_buckets[MAX_BUCKET_NUM];

  Direction *direction;


 private:
  Molecule *m_molecule;

  double m_radius;
  double m_bucketSize;
  int m_numBuckets; //number of buckets

  int m_numDOFs;
  std::list<Configuration *> m_samples;
  std::vector<Configuration *> m_path;

  double m_top_min_rmsd;
  int m_top_min_rmsd_id;
  double m_minMovDihDistance;
  int m_minMovDihDistance_id;

  int m_numSamples;
  int m_gradientSelection;
  double m_stepSize;
  double m_maxRotation;
  bool m_scaleToRadius;
};


#endif

