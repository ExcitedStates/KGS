/*KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data*/

/* DihedralRRT.cpp
 *
 *  Created on: 19.08.2015
 *      Author: StDoBudd
 */

#include "DihedralRRT.h"

#include <iomanip>
#include <stack>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "Logger.h"
#include "core/Transformation.h"
#include "CTKTimer.h"
#include "metrics/RMSD.h"
#include "metrics/Dihedral.h"


using namespace std;

//const double RMSD_DIFF_THRESHOLD = 0.01;
const double MOV_DIH_THRESHOLD = 0.001;

extern double selectNodeTime;



DihedralRRT::DihedralRRT(Molecule *protein, Move& move, Direction& direction):
	SamplingPlanner(move),
	m_max_distance(SamplingOptions::getOptions()->explorationRadius),
  direction(direction)
{
	m_protein = protein;
	m_numDOFs = m_protein->m_spanning_tree->getNumDOFs();
	Configuration *pSmp = new Configuration(m_protein);
	pSmp->updateProtein();
  pSmp->computeCycleJacobianAndNullSpace();
	m_protein->m_conf = pSmp;
	m_protein->m_conf_backup = pSmp;
	m_target = NULL;
	m_samples.push_back(pSmp);
	pSmp->m_vdwEnergy = 99999;
	pSmp->m_id = 0; // root
	if(SamplingOptions::getOptions()->metric_string=="rmsd"){
		cerr<<"Metric option rmsd incompatible with planning strategy dihedralRRT. Please use dihedral m_metric."<<endl;
		exit(-1);
	}
	if(SamplingOptions::getOptions()->metric_string=="dihedral") 	metric = new metrics::Dihedral();

	//RMSD for observation purposes
	rmsdMetric = new metrics::RMSD();

	m_deform_mag = 0.25;
	m_rand_radius = 2;

	m_top_min_rmsd = 99999;
	m_top_min_rmsd_id = -1;
	m_numSamples = 0;
	m_minMovDihDistance = 99999;
	m_minMovDihDistance_id = -1;

}

DihedralRRT::~DihedralRRT() {
	Configuration *pSmp;
	for (list<Configuration *>::iterator iter = m_samples.begin(); iter != m_samples.end(); iter++)
	{
		delete *iter;
	}
	
}



void DihedralRRT::GenerateSamples(){
	int nBatch = SamplingOptions::getOptions()->samplesToGenerate;
	int sample_id=0, max_depth=0, failed_trials=0, total_trials=0;
	Configuration *pTarget=NULL, *pClosestSmp, *pNewSmp=NULL;
  gsl_vector* gradient = gsl_vector_alloc(m_protein->totalDofNum());

	//Save initial file (for movie)
	writeNewSample(m_samples.front(), m_samples.front(), sample_id);

	bool createNewTarget = false;

	while (sample_id<nBatch) {
		++total_trials;

		CTKTimer timer;
		double start_time = timer.getTimeNow();

		if( SamplingOptions::getOptions()->sampleRandom || pNewSmp == NULL || createNewTarget) {
			log("dominik")<<"Generating new target, getting new seed"<<endl;
			pTarget = GenerateRandConf(); // used in selection ONLY if no m_target m_protein specified
			createNewTarget = false;
			pClosestSmp = SelectNodeFromBuckets(pTarget);
			double end_time = timer.getTimeNow();
			selectNodeTime += end_time - start_time;
		}else{
			log("dominik")<<"Using latest sample as seed"<<endl;
			pClosestSmp = m_samples.back();
		}

    direction.gradient(pClosestSmp, NULL, gradient);
		pNewSmp = move.move(pClosestSmp, gradient);

		if( !pNewSmp->updatedProtein()->inCollision() ) {
			++sample_id;
			m_numSamples=sample_id;

			pNewSmp->m_distanceToIni = metric->distance(pNewSmp,m_samples.front());
			pNewSmp->m_distanceToParent = metric->distance(pNewSmp,pClosestSmp);
      pNewSmp->m_id = sample_id;

			writeNewSample(pNewSmp, m_samples.front(), sample_id);

			if (pNewSmp->m_treeDepth > max_depth)
				max_depth = pNewSmp->m_treeDepth;

			double distToRandGoal = metric->distance(pNewSmp,pTarget);

			//log("samplingStatus") << "> New structure: newpdb_"<<sample_id<<".pdb .. RMSD to initial: "<< pNewSmp->m_rmsd_initial<<endl;
			log("samplingStatus") << "> New structure: " << SamplingOptions::getOptions()->moleculeName<<"_new_" << sample_id << ".pdb";
			log("samplingStatus") << " .. Distance to initial: " << setprecision(6) << pNewSmp->m_distanceToIni;
			log("samplingStatus") << " .. Distance to current target: " << setprecision(3) << distToRandGoal;
			log("samplingStatus") << " .. Null-space dimension: " << pNewSmp->getNullspace()->NullspaceSize()<<endl;

			if(distToRandGoal <= MOV_DIH_THRESHOLD){//current target reached
				delete pTarget;
				createNewTarget = true;
			}
			//Check if target has been reached!
			//if(m_target != NULL && (m_top_min_rmsd < RMSD_DIFF_THRESHOLD )){
			//cout<<"Reached target, not creating more samples!"<<" rmsd: "<<m_top_min_rmsd<<", dih: "<<m_minMovDihDistance<<endl;
			//break;
			//}
		} else {
			++failed_trials;
			//Only for the straight path testing
			//if( !options.sampleRandom )
			//return NULL;
			delete pTarget;
			createNewTarget = true;
		}
	}

  gsl_vector_free(gradient);
}

Configuration* DihedralRRT::GenerateRandConf() {
	//Configuration *pNewSmp = new Configuration(getNumDOFs());
	Configuration *pNewSmp = new Configuration(m_protein);

	double length=0.0;
	for (int i=0; i<m_numDOFs; ++i) {
		pNewSmp->m_dofs[i] = dPi * RandomN1P1();
		length+= pNewSmp->m_dofs[i] * pNewSmp->m_dofs[i];
	}
	length = sqrt(length);
	if(SamplingOptions::getOptions()->scaleToRadius){
		double factor = pow( Random01(), 1.0/m_numDOFs )*SamplingOptions::getOptions()->explorationRadius/length;
		for (int i=0; i<m_numDOFs; ++i) {
			pNewSmp->m_dofs[i] = factor * pNewSmp->m_dofs[i];
		}
	}
	
	pNewSmp->m_id = -1;
	m_protein->SetConfiguration(pNewSmp);
	pNewSmp->updateGlobalTorsions();
	return pNewSmp;
}

Configuration *DihedralRRT::SelectNodeFromBuckets (Configuration *pTarget) {
	Configuration *pSmp, *pMinSmp;
	double min_distance = 1000000.0;
	double distance;
	pMinSmp = NULL;

	for (list<Configuration*>::iterator iter=m_samples.begin(); iter!=m_samples.end(); ++iter) {
		pSmp = *iter;
		distance = metric->distance(pSmp,pTarget);
		if (distance < min_distance) {
			min_distance = distance;
			pMinSmp = pSmp;
		}
	}

	return pMinSmp;
}



