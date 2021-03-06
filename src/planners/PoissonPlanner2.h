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


#ifndef POISSONPLANNER2_H_
#define POISSONPLANNER2_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "metrics/Metric.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"

class Molecule;

typedef std::tuple<Residue *, Residue *, Residue *> ResTriple;

/**
 * A sampling planner based on Poisson-disc sampling as described in Robert Bridsons 2007 SIGGRAPH paper.
 * Each new sample is generated by randomly picking an existing open 'seed' sample, make a predefined number of
 * attempts to generate valid perturbations, adde them to the open set, and finally close the seed sample. A
 * valid perturbation gives a non-clashing sample that is not too close to existing samples.
 *
 * This is repeated until stop_after samples have been generated, or until there are no more 'open' samples.
 */
class PoissonPlanner2 : public SamplingPlanner {
 public:
  PoissonPlanner2(Molecule *, std::vector<ResTriple> &, int, int, double, const std::string &);

  ~PoissonPlanner2();

  void generateSamples();

  std::list<Configuration *> &getSamples() { return all_samples; }

  bool m_checkAll = false;
 private:
  /// Even if open_samples is non-empty sampling will stop after this many new samples have been generated
  const int m_stopAfter;
  const int m_maxRejectsBeforeClose; ///< Number perturbations that are tried before a sample is 'closed'
  const double m_bigRad;              ///< Largest allowed step-size
  const double m_lilRad;              ///< Smallest allowed distance between any two samples
  const std::vector<ResTriple> &m_ikTriples; ///< Triples of residues that will be rebuilt before closing.
  const std::string m_resNetwork;

  std::list<Configuration *> open_samples;     ///< Non-closed samples
  std::list<Configuration *> closed_samples;   ///< Samples that have tested more than max_rejects_before_close perturbations
  std::list<Configuration *> all_samples;      ///< For convenience and return

  std::vector<std::tuple<Residue *> > m_tripeptides; ///< Preprocessed residue triples for use in exact IK

  Molecule *m_protein;

  Configuration *m_root;

  /** The largest distance from a sample to any of its descendants */
  std::map<Configuration *, double> m_maxDist;

  /** Memoized distances between configuration pointers */
  std::map<std::pair<Configuration *, Configuration *>, double> m_distances;

  /** Compute distance between c1 and c2 using the memoized map if possible */
  double memo_distance(Configuration *c1, Configuration *c2);

  /** Collect all open or closed configurations within dist of conf and add them to ret. */
  void collectPossibleChildCollisions(Configuration *conf,
                                      std::vector<Configuration *> &ret,
                                      double childOffset);

  void collectPossibleChildCollisions(Configuration *conf,
                                      std::vector<Configuration *> &ret,
                                      Configuration *v,
                                      double childOffset);

  void updateMaxDists(Configuration *newConf);

  void updateMaxDists(Configuration *v, Configuration *newConf);
};

#endif /* POISSONPLANNER_H_ */
