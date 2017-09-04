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

#ifndef MCMCPLANNER_H_
#define MCMCPLANNER_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <directions/Direction.h>

#include "metrics/Metric.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"

class Molecule;

/**
 * A sampling planner based on MCMC sampling.
 */
class MCMCPlanner : public SamplingPlanner {
 public:
  MCMCPlanner(
      Molecule *molecule,
      Direction *direction,
      int stopAfter,
      double stepSize
  );

  ~MCMCPlanner();

  void generateSamples();

  std::list<Configuration *> &getSamples() { return m_samples; }

 private:
  std::list<Configuration *> m_samples;      ///< For convenience and return

  const int m_stopAfter;              ///< Stop after this many samples have been generated
  const double m_bigRad;              ///< Largest allowed step-size
  const double m_lilRad;              ///< Smallest allowed distance between any two samples
  Direction* m_direction;




  Molecule *m_molecule;
};

#endif /* MCMCPLANNER_H_ */
