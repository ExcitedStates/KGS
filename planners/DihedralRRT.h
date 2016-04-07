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

#include "directions/Direction.h"
#include "metrics/Metric.h"
#include "SamplingOptions.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"


class Molecule;
class Configuration;

#define DEFAULT_MAX_RMSD 25



class DihedralRRT : public SamplingPlanner{
public:
	DihedralRRT(Molecule *, Move&, Direction&);
	~DihedralRRT();

	void GenerateSamples();

	ConfigurationList& Samples(){ return m_samples; }

	double m_deform_mag;
	double m_rand_radius;

protected:
	Configuration* GenerateRandConf();
	Configuration* SelectNodeFromBuckets(Configuration *pTarget);

	metrics::Metric* metric;
	metrics::Metric* rmsdMetric;

  Direction& direction;


public:
	Molecule *m_protein;
	Molecule *m_target;

	double m_max_distance;

	ConfigurationList m_samples;
	int m_DOF;
	ConfigurationArray m_path;

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
