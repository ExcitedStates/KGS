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
#ifndef RRTPLANNER_H
#define RRTPLANNER_H

#include <iostream>
#include <fstream>
#include <string>
#include <directions/Direction.h>

#include "metrics/Metric.h"
#include "SamplingOptions.h"
#include "core/Configuration.h"
#include "planners/SamplingPlanner.h"


class Molecule;

#define MAX_BUCKET_NUM 101
#define NUM_BINS 100
#define DEFAULT_MAX_RMSD 25


class RRTPlanner: public SamplingPlanner
{
public: 
	RRTPlanner(Molecule *, Move& move, metrics::Metric& metric, Direction& direction);
	~RRTPlanner();

	void GenerateSamples();

	ConfigurationList& Samples(){ return m_samples; }

	double m_deform_mag;
	double m_rand_radius;

protected:
	Configuration* GenerateRandConf();
	Configuration* SelectNodeFromBuckets(Configuration *pTarget);

	unsigned int m_current_max_bucket_id;
	std::list<Configuration*> distance_buckets[MAX_BUCKET_NUM];

  Direction& direction;


public:
	Molecule *m_protein;
	Molecule *m_target;

	double m_max_distance;
	double m_bucket_size;
	int m_numBuckets; //number of buckets

	int m_numDOFs;
	ConfigurationList m_samples;
	ConfigurationArray m_path;

	double m_top_min_rmsd;
	int m_top_min_rmsd_id;
	double m_minMovDihDistance;
	int m_minMovDihDistance_id;

	int m_numSamples;
};


#endif

