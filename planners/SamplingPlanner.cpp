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

  std::list<Configuration *> &m_samples = Samples();

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

