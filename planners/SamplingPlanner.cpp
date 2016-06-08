#include <math/gsl_helpers.h>

#include "SamplingPlanner.h"

#include "IO.h"

using namespace std;

double selectNodeTime = 0;

SamplingPlanner::SamplingPlanner(Move& move, metrics::Metric& metric):
	move(move),
	m_metric(metric)
{

}

SamplingPlanner::~SamplingPlanner(){}

void SamplingPlanner::createTrajectory(){

	std::list<Configuration*>::iterator cit;
	Configuration *pSmp;

	std::list<Configuration*>& m_samples = Samples();

//	if( !( (SamplingOptions::getOptions()->targetStructureFile).empty() ) ){ //we have a target, trajectory to the closest configuration
//		for( cit=m_samples.begin(); cit!=m_samples.end(); cit++ ){
//			if((*cit)->m_id == m_closestFwdSample->m_id ){
//				pSmp = (*cit);
//				break;
//			}
//		}
//	}
//	else{//no target, use longest path
		int treeDepth = 0;
		for( cit=m_samples.begin(); cit!=m_samples.end(); cit++ ){
				if((*cit)->m_treeDepth >= treeDepth){//m_top_min_rmsd_id
					treeDepth = (*cit)->m_treeDepth;
					pSmp = (*cit);
				}
			}
	//}
	//m_molecule->setConfiguration(pSmp);
	Molecule * m_protein = pSmp->updatedMolecule();
  SamplingOptions& options = *(SamplingOptions::getOptions());

	const string& out_path = options.workingDirectory;
	const string& name = m_protein->getName();

	string out_collPdb = out_path + "output/" + name + "_path.pdb";
	///save pyMol movie script
	string out_pyMol=out_path + "output/" +  name + "_pathMov.pml";

	IO::writeTrajectory(m_protein, out_collPdb, out_pyMol);
}

void SamplingPlanner::writeNewSample(Configuration* conf, Configuration* ref, int sample_num)
{
	const string& out_path = SamplingOptions::getOptions()->workingDirectory;
	const string& name = conf->getMolecule()->getName();
	string out_file = out_path+"output/"+name+"_new_"+std::to_string(static_cast<long long>(sample_num))+".pdb";

	if(SamplingOptions::getOptions()->saveData > 0){

		Molecule * protein = conf->updatedMolecule();
		IO::writePdb(protein, out_file);
	}


	if(SamplingOptions::getOptions()->saveData > 1){
		Molecule * protein = conf->updatedMolecule();

		string out_q = out_path+"output/"+name+"_q_"+std::to_string(static_cast<long long>(sample_num))+".txt";

		IO::writeQ(protein, ref, out_q);
	}


	if(SamplingOptions::getOptions()->saveData > 3){
		Molecule * protein = conf->updatedMolecule();

		// Save Jacobian and Nullspace to file
		string outJac=out_path + "output/" +  name + "_jac_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";
		string outNull=out_path + "output/" +  name + "_nullSpace_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";
		// Save singular values
		string outSing=out_path + "output/" +  name + "_singVals_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";
		// Save pyMol coloring script
		string pyMol=out_path + "output/" +  name + "_pyMol_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".pml";
		string rbFile=out_path + "output/" +  name + "_RBs_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";
		string covFile=out_path + "output/" +  name + "_covBonds_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";
		string statFile=out_path + "output/" +  name + "_stats_" +
			std::to_string(static_cast<long long>(sample_num))
			+ ".txt";

		IO::writePyMolScript(protein, out_file, pyMol);

		///Write Jacobian
		//gsl_matrix_outtofile(conf->CycleJacobian, outJac);
		///Write Null Space matrix
		//gsl_matrix_outtofile(conf->CycleNullSpace->V,outNull);
		//gsl_vector_outtofile(conf->CycleNullSpace->singularValues,outSing);

		conf->getNullspace()->WriteMatricesToFiles(outJac, outNull, outSing);
		IO::writeRBs(protein, rbFile);
		IO::writeStats(protein,statFile);

	}

}
