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
#ifndef SAMPLINGPLANNER_H_
#define SAMPLINGPLANNER_H_

#include "core/Configuration.h"
#include "moves/Move.h"
#include "metrics/Metric.h"

/**
 * The superclass for all sampling planners.
 * A sampling planner encodes a strategy for generating a list of samples.
 *
 * Each returned configuration should have a unique id (Configuration::m_ID).
 */
class SamplingPlanner{
 public:
  //SamplingPlanner(Move& move, metrics::Metric& metric);
  SamplingPlanner();

  virtual ~SamplingPlanner() = 0;

  virtual void initialize(Move* move, metrics::Metric* metric, const std::string& workingDir, int saveData);

  /** Generate samples. */
  virtual void generateSamples() = 0;

  /** Return a reference to the list of all generated samples */
  virtual std::list<Configuration*>& getSamples() = 0;
  
  /**
   * Depending on the options, either the sample at the end of the longest path
   * or the sample nearest to the goal conformation is chosen, and the path from
   * initial to this sample is written to a file.
   */
  virtual void createTrajectory();

//  void writeNewSample(Configuration* conf, Configuration* ref, int sample_num);

 protected:

  Move* m_move;                ///< Performs conformational perturbations.
  metrics::Metric* m_metric;   ///< Metric used to measure distances between samples
  std::string m_workingDir;    ///< Working directory into which samples should be written
  int m_saveData;              ///< Amount of data that should be written during sampling



};

#endif

