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

#include <math/gsl_helpers.h>
#include <math/NullspaceSVD.h>

#include "SamplingPlanner.h"

#include "IO.h"

using namespace std;

double selectNodeTime = 0;

SamplingPlanner::SamplingPlanner() :
    m_move(nullptr),
    m_metric(nullptr),
    m_workingDir(""),
    m_saveData(-1) {
}

SamplingPlanner::~SamplingPlanner() {}

void SamplingPlanner::initialize(Move *move, metrics::Metric *metric, const std::string &workingDir, int saveData) {
  m_move = move;
  m_metric = metric;
  m_workingDir = workingDir;
  m_saveData = saveData;
}

void SamplingPlanner::createTrajectory() {
  if (m_move == nullptr) {
    cerr << "SamplingPlanner::createTrajectory - ERROR: SamplingPlanner must be initialized before use" << endl;
    exit(-1);
  }

  std::list<Configuration *>::iterator cit;
  Configuration *pSmp;

  std::list<Configuration *> &m_samples = getSamples();

//	if( !( (ExploreOptions::getOptions()->targetStructureFile).empty() ) ){ //we have a target, trajectory to the closest configuration
//		for( cit=m_samples.begin(); cit!=m_samples.end(); cit++ ){
//			if((*cit)->m_id == m_closestFwdSample->m_id ){
//				pSmp = (*cit);
//				break;
//			}
//		}
//	}
//	else{//no target, use longest path
  int treeDepth = 0;
  for (cit = m_samples.begin(); cit != m_samples.end(); cit++) {
    if ((*cit)->m_treeDepth >= treeDepth) {//m_top_min_rmsd_id
      treeDepth = (*cit)->m_treeDepth;
      pSmp = (*cit);
    }
  }
  //}
  //m_molecule->setConfiguration(pSmp);
  Molecule *m_protein = pSmp->updatedMolecule();

//	ExploreOptions& options = *(ExploreOptions::getOptions());
//	const string& out_path = options.workingDirectory;
  const string &out_path = m_workingDir;
  const string &name = m_protein->getName();

  string out_collPdb = out_path + "output/" + name + "_path.pdb";
  ///save pyMol movie script
  string out_pyMol = out_path + "output/" + name + "_pathMov.pml";

  IO::writeTrajectory(m_protein, out_collPdb, out_pyMol);
}

