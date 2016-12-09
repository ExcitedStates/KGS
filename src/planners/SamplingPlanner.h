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

