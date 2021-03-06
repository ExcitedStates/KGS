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



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "metrics/RMSD.h"
#include "IO.h"
#include "core/Chain.h"

using namespace std;

Molecule * readProtein(char* path){
	char* tmp = realpath(path, nullptr);
	if(tmp==nullptr){ cerr<<path<<" is not a valid PDB-file"<<endl; exit(-1); }

	Selection movingResidues("all");
	Molecule* protein = IO::readPdb(path);
	protein->initializeTree(movingResidues,1.0);
	return protein;
}

double torsion(Bond * bond){
	Atom *a1,*a2,*a3,*a4;
	a2 = bond->m_atom1;
	a3 = bond->m_atom2;
	//Find min neighbor of a2 (not a3)
	a1 = a2->Cov_neighbor_list[0]; if(a1==a3) a1 = a2->Cov_neighbor_list[1];
	for(int i=1;i<a2->Cov_neighbor_list.size();i++) if( a2->Cov_neighbor_list[i]!=a3 && a2->Cov_neighbor_list[i]->getName()<a1->getName() ) a1 = a2->Cov_neighbor_list[i];
	a4 = a3->Cov_neighbor_list[0]; if(a4==a2) a4 = a3->Cov_neighbor_list[1];
	for(int i=1;i<a3->Cov_neighbor_list.size();i++) if( a3->Cov_neighbor_list[i]!=a2 && a3->Cov_neighbor_list[i]->getName()<a4->getName() ) a4 = a3->Cov_neighbor_list[i];

	double torsion = TorsionalAngle(a1->m_position, a2->m_position, a3->m_position, a4->m_position);
	return torsion;
}

void collectConfigurations(Molecule * native, int arrSz, char* fileList[], vector<Configuration*> &confs, vector<Molecule *> &prots){

	for(int i=0;i<arrSz;i++){
		Molecule * struc = readProtein(fileList[i]);
		Configuration* conf = new Configuration(native);
		//Configuration* conf = new Configuration(native->m_spanningTree->getNumDOFs());

    /*
		for(auto vit = struc->m_spanningTree->m_vertices.begin(); vit != struc->m_spanningTree->m_vertices.end(); vit++){
			KinVertex* vertex = vit->second;
			if(vertex->isRibose) {
				SugarVertex* v = reinterpret_cast<SugarVertex*>(vertex);
				double strucTorsion = v->initTorsion;
				double nativeTorsion = 1000;

				for(auto nvit = native->m_spanningTree->m_vertices.begin(); nvit != native->m_spanningTree->m_vertices.end(); nvit++){
					KinVertex* nvertex = nvit->second;
					if(nvertex->isRibose) {
						SugarVertex* nv = reinterpret_cast<SugarVertex*>(nvertex);
						if(nv->getDOF()->getIndex()==v->getDOF()->getIndex()){
							nativeTorsion=nv->initTorsion;
							break;
						}
					}
				}

				double diff = strucTorsion-nativeTorsion;
				if(diff<-CTK_PI) diff+=2*CTK_PI;
				if(diff>CTK_PI) diff-=2*CTK_PI;
				//cout<<diff<<" ";
				if(diff>3.15){
					cout<<"Sugar: ";
					cout<<v->getDOF()->getIndex()<<" .. "<<diff<<" .. "<<strucTorsion<<" .. "<<nativeTorsion<<endl;
					exit(-1);
				}
				conf->m_dofs[v->getDOF()->getIndex()] = diff;

			}
		}
     */

		for(vector<KinEdge*>::iterator eit=struc->m_spanningTree->m_edges.begin(); eit!=struc->m_spanningTree->m_edges.end(); ++eit){
			KinEdge* e = *eit;

			double strucTorsion = torsion(e->getBond());
			double nativeTorsion = 1000;

			for(vector<KinEdge*>::iterator neit=native->m_spanningTree->m_edges.begin(); neit!=native->m_spanningTree->m_edges.end(); ++neit){
				KinEdge* ne = *neit;
				if(ne->getDOF()->getIndex()==e->getDOF()->getIndex()){
					nativeTorsion = torsion(ne->getBond());
					break;
				}
			}
			double diff = strucTorsion-nativeTorsion;
			if(diff<-CTK_PI) diff+=2*CTK_PI;
			if(diff>CTK_PI) diff-=2*CTK_PI;
			//cout<<diff<<" ";
			if(diff>3.15){
					cout<<"KinEdge: ";
				cout<<e->getDOF()->getIndex()<<" .. "<<diff<<" .. "<<strucTorsion<<" .. "<<nativeTorsion<<endl;
				exit(-1);
			}
			conf->m_dofs[e->getDOF()->getIndex()] = diff;

		}

		confs.push_back(conf);
		prots.push_back(struc);
		cout<<".";cout.flush();

	}
}

int main(int argc, char* argv[]){
	if(argc<3){
        cout<<"Usage: "<<argv[0]<<" <native pdb> <pdb-list>"<<endl;
        cout<<"Computes PCA of degrees of freedom. Outputs deformability and mobility as b-factors."<<endl;
        return -1;
    }
	cout<<"Using "<<argv[1]<<" as native (wont be included in ensemble PCA)."<<endl;

	cout<<"Collecting configurations "; cout.flush();
	Molecule * native = readProtein(argv[1]);
	//native->BackupAtomPos();

	vector<Configuration*> configurations;
	vector<Molecule *> proteins;
	collectConfigurations(native, argc-2, &(argv[2]), configurations, proteins);
	cout<<"done. Total: "<<configurations.size()<<endl;

	cout<<"Calculating covariance matrix .. ";
	int n = native->m_spanningTree->getNumDOFs();
	double* cov = new double[n*n];
	double* avg = new double[n];
	for(int i=0;i<n;i++){ avg[i]=0; for(int j=0;j<n;j++) cov[i*n+j] = 0; }

	for(vector<Configuration*>::iterator cit=configurations.begin(); cit!=configurations.end(); cit++){
		for(int i=0;i<n;i++){
			double val = (*cit)->m_dofs[i];
			if(val>3.1416){
				cerr<<"Error: dof["<<i<<"] = "<<val<<endl;
				exit(-1);
			}
			avg[i]+= val;
		}
	}
	for(int i=0;i<n;i++){ avg[i]/=configurations.size();}
	for(int i=0;i<n;i++)
		for(int j=i;j<n;j++){
			for(vector<Configuration*>::iterator cit=configurations.begin(); cit!=configurations.end(); cit++){
				cov[i*n+j] += ((*cit)->m_dofs[i] - avg[i]) * ((*cit)->m_dofs[j] - avg[j]);
			}
			cov[i*n+j]/=configurations.size();
			cov[j*n+i]=cov[i*n+j];
		}
	cout<<" done"<<endl;


	cout<<"Calculating eigenvectors .. ";cout.flush();
	gsl_matrix_view m = gsl_matrix_view_array (cov, n, n);

	gsl_vector *eval = gsl_vector_alloc (n);
	gsl_matrix *evec = gsl_matrix_alloc (n, n);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
	gsl_eigen_symmv (&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	cout<<"done"<<endl;


		for (int i = 0; i < n; i++)
		{
			double eval_i = gsl_vector_get (eval, i);
			printf ("eigenvalue = %g\n", eval_i);

			//gsl_vector_view evec_i = gsl_matrix_column (evec, i);
			//printf ("eigenvector = \n");
			//gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
		}

	double atom_def[native->getAtoms().size()];
	double atom_mob[native->getAtoms().size()];
	for(int i=0;i<native->getAtoms().size(); i++){
		atom_def[i] = 0;
		atom_mob[i] = 0;
	}

	for(int component=0;component<evec->size2;component++){
		gsl_vector_view evec_i = gsl_matrix_column (evec, component);
		double eval_i = gsl_vector_get(eval, component);
		
		for(vector<KinEdge*>::iterator eit = native->m_spanningTree->m_edges.begin(); eit != native->m_spanningTree->m_edges.end(); ++eit){
			KinEdge* e = *eit;
			int dof = e->getDOF()->getIndex();
			double evec_component = gsl_vector_get(&(evec_i.vector), dof);
			evec_component = fabs(evec_component)*eval_i*eval_i*100.0;
			int id1 = e->getBond()->m_atom1->getId();
			int id2 = e->getBond()->m_atom2->getId();
			atom_def[id1]+=evec_component;
			atom_def[id2]+=evec_component;
//			cout<<e<<" "<<evec_component<<endl;
		}

    /*
		for (auto vit=native->m_spanningTree->m_vertices.begin(); vit!=native->m_spanningTree->m_vertices.end(); vit++){
			if( (*vit).second->isRibose ){
				SugarVertex* v = reinterpret_cast<SugarVertex*>((*vit).second);
				int dof = v->getDOF()->getIndex();
				double evec_component = gsl_vector_get(&(evec_i.vector), dof);
				evec_component = fabs(evec_component)*eval_i*eval_i;
				//cout<<v<<" "<<evec_component<<endl;
				for(auto ait=v->m_rigidbody->Atoms.begin(); ait!=v->m_rigidbody->Atoms.end(); ++ait){
					int aid = (*ait)->getId();
					//atom_def[aid]+=evec_component*0.22*eval_i*eval_i*100.0;
					atom_def[aid]+=evec_component*100.0;

				}

			}
		}
     */
	}

//    for(int i=0;i<configurations.size(); i++){
//        Molecule* prot = proteins[i];
//        RRTPlanner::align(prot, native, nullptr);
//        for (int a=0;a<native->atoms.size();a++){
//            Atom* natom = native->atoms[a];
//            Atom* atom = prot->atoms[a];
//            atom_mob[a]+= atom->m_position.distance(natom->m_position)/configurations.size();
//		}
//	}
    double defMax = 0.0;
//    double mobMax = 0.0;
	for(int i=0;i<native->getAtoms().size(); i++){
        if(atom_def[i]>defMax) defMax = atom_def[i];
//        if(atom_mob[i]>mobMax) mobMax = atom_mob[i];
	}
	cout<<"DefMax: "<<defMax<<endl;
	for(int i=0;i<native->getAtoms().size(); i++){
        atom_def[i] /= defMax;
//        atom_mob[i] /= mobMax;
	}

	ofstream output("pca_def.pdb");

  for(auto const& atom: native->getAtoms()){
    Residue* res = atom->getResidue();
    char buffer[100];
    sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ",
            atom->getId(),atom->getName().c_str(),
            res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
            atom->m_position.x,atom->m_position.y,atom->m_position.z,atom_def[atom->getId()],atom_def[atom->getId()], atom->getType().c_str());
    string line(buffer);
    output << line << endl;
  }
  output.close();

//	ofstream output2("pca_mob.pdb");
//	//for (vector<Atom*>::iterator atom_itr=native->atoms.begin(); atom_itr!=native->atoms.end(); ++atom_itr) {
//    for (int a=0;a<native->atoms.size(); a++){
//        Atom* atom = native->atoms[a]; // *atom_itr;
//		Residue* res = atom->getResidue();
//		char buffer[100];
//		sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  ",
//			atom->getId(),atom->Original_name.c_str(),
//			res->getName().c_str(),res->getChainId().c_str(),res->getId(),
//			atom->m_position.x,atom->m_position.y,atom->m_position.z,atom_mob[atom->getId()],atom_mob[atom->getId()], atom->getType().c_str());
//		string line(buffer);
//		output2 << line << endl;
//	}
//	output2.close();

    double amplitude = 1.0;
    double stepsize = 0.2;

	if(true) return 0;
	
	for(int component=0;component<3;component++){
		cout<<"Writing component "<<component<<" pdb files ";cout.flush();
		Molecule * aligned = readProtein(argv[1]);
		Configuration* conf = new Configuration(aligned);
		//Configuration* conf = new Configuration(native->m_spanningTree->getNumDOFs());
		gsl_vector_view evec_i = gsl_matrix_column (evec, component);
		for(double a=-amplitude;a<=amplitude;a+=stepsize){
            if(fabs(a)<0.001) a = 0.0;
			for(int i=0;i<n;i++) {
				conf->m_dofs[i] = avg[i] + gsl_vector_get(&(evec_i.vector), i) * a;
			}
      aligned->setConfiguration(conf);
			//RRTPlanner::align(aligned, native, nullptr);
      Selection sel;

			metrics::RMSD rmsd(sel);
      rmsd.align(aligned, native);
			stringstream ss;ss<<"comp_"<<component<<"_"<<a<<".pdb";
			IO::writePdb(aligned, ss.str());
			cout<<".";cout.flush();
		}
		delete conf;
		cout<<" done"<<endl;


		stringstream combineCommand;
		combineCommand<<"trunk/Scripts/combinePDBs.py ";
		for(double a=0;a<amplitude;a+=stepsize){
            if(fabs(a)<0.001) a = 0.0;
			stringstream ss;ss<<"comp_"<<component<<"_"<<a<<".pdb";
			combineCommand<<ss.str()<<" ";
		}
		for(double a=amplitude;a>-amplitude;a-=stepsize){
            if(fabs(a)<0.001) a = 0.0;
			stringstream ss;ss<<"comp_"<<component<<"_"<<a<<".pdb";
			combineCommand<<ss.str()<<" ";
		}
		for(double a=-amplitude;a<0;a+=stepsize){
            if(fabs(a)<0.001) a = 0.0;
			stringstream ss;ss<<"comp_"<<component<<"_"<<a<<".pdb";
			combineCommand<<ss.str()<<" ";
		}
		combineCommand<<" > comp_"<<component<<".pdb";
		cout<<"To create multi-model PDB, run:"<<endl;
		cout<<combineCommand.str()<<endl;
	}
	gsl_vector_free (eval);
	gsl_matrix_free (evec);

}
