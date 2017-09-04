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

#include "MCMCPlanner.h"

#include <iomanip>
#include <stack>
#include <IO.h>

#include "Logger.h"
#include "core/Transformation.h"

using namespace std;

MCMCPlanner::MCMCPlanner(
    Molecule * molecule,
    Direction * direction,
    int stopAfter,
    double stepSize
):
  SamplingPlanner(),
  m_molecule(molecule),
  m_direction(direction),
  m_stopAfter(stopAfter),
  m_bigRad(stepSize*4.0/3.0),
  m_lilRad(m_bigRad/2)
{
}

MCMCPlanner::~MCMCPlanner() {
  for (auto &s : m_samples) {
    delete s;
  }
}

void MCMCPlanner::generateSamples()
{
  log("samplingStatus")<<"Sampling ...\n";
  int clashes = 0;

  gsl_vector* gradient = gsl_vector_alloc(m_molecule->m_spanningTree->getNumDOFs());
  m_samples.push_back(new Configuration(m_molecule));
  for (int i = 0; i < m_stopAfter; i++) {
    Configuration* seed = m_samples.back();
    m_direction->gradient(seed, nullptr, gradient);
    //gsl_vector_scale_max_component(gradient, options.maxRotation);

    // Scale gradient so move is in Poisson disc
    Configuration *new_conf = m_move->move(seed, gradient); //Perform move
    double dist = m_metric->distance(new_conf, seed);
    int scaleAttempts = 0;
    while( dist<m_lilRad || dist>m_bigRad){
      if(++scaleAttempts==5) break;
      double gradientScale = (m_bigRad+m_lilRad)/(2.0*dist);
      gsl_vector_scale(gradient, gradientScale);
      delete new_conf;
      new_conf = m_move->move(seed, gradient);
      dist = m_metric->distance(new_conf, seed);
    }

//      scale_gradient(gradient, &protein);
//      gsl_vector_scale(gradient, options.stepSize);
//      Configuration* new_conf = move->move(seed, gradient);

    if(new_conf->updatedMolecule()->inCollision()) {
      clashes++;
      i--;
    } else {
      log("samplingStatus")<<"Accepted conformation "<<i<<endl;
      IO::writePdb(new_conf->updatedMolecule(), "output/conf_" + std::to_string((long long) i) + ".pdb");
      m_samples.push_back(new_conf);
    }
  }


  //Print final status
  log("samplingStatus")<<"MCMC-planner: Rejects from clash:          "<<clashes<<endl;

  log("samplingStatus")<<"Done"<<endl;
}


