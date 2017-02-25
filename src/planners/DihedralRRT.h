/*
 * DihedralRRT.h
 *
 *  Created on: 19.08.2015
 *      Author: StDoBudd
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
      bool sampleRandom
  );

  ~DihedralRRT();

  void generateSamples();

  std::list<Configuration *> &getSamples() { return m_samples; }

  double m_deform_mag;
  double m_rand_radius;

 protected:
  Configuration *GenerateRandConf();

  Configuration *SelectClosestNode(Configuration *pTarget);

  Direction *m_direction;

 private:
  Molecule *m_protein;
  Molecule *m_target;

  metrics::Dihedral *m_dihedralMetric;

  double m_maxDistance;
  bool m_sampleRandom;

  std::list<Configuration *> m_samples;
  int m_numDOFs;
  std::vector<Configuration *> m_path;

  int m_nCDCall;
  int m_nRMSDCall;

  double m_top_min_rmsd;
  int m_top_min_rmsd_id;
  double m_minMovDihDistance;
  int m_minMovDihDistance_id;

  int m_numSamples;
  int m_max_depth;
};

#endif /* DIHEDRALRRT_H_ */
