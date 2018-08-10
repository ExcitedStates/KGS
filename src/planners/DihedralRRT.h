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


#ifndef DIHEDRALRRT_H_
#define DIHEDRALRRT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <metrics/Dihedral.h>

#include "directions/Direction.h"
#include "metrics/Metric.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"


#define DEFAULT_MAX_RMSD 25


class DihedralRRT : public SamplingPlanner {
 public:
  DihedralRRT(
      Molecule *protein,
      Direction *direction,
      int numSamples,
      double maxDistance,
      double maxRotation,
      bool sampleRandom
  );

  ~DihedralRRT();

  void generateSamples();

  std::list<Configuration *> &getSamples() { return m_samples; }

  Configuration* getFurthestSample() {return m_furthestSample;}

 protected:
  Configuration *GenerateRandConf();

  Configuration *SelectClosestNode(Configuration *pTarget);

  Direction *m_direction;

 private:
  Molecule *m_protein;
  Molecule *m_target;

  metrics::Dihedral *m_dihedralMetric;

  double m_explorationRadius;
  bool m_sampleRandom;

  std::list<Configuration *> m_samples;
  int m_numDOFs;
  std::vector<Configuration *> m_path;

  int m_numSamples;
  double m_maxRotation;
  int m_max_depth;

  Configuration* m_furthestSample;
  double m_furthestDistance;
};

#endif /* DIHEDRALRRT_H_ */
