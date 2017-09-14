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

#include "TestDummy.h"

#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h> //system

#include "../Logger.h"
#include "core/Molecule.h"
#include "../IO.h"
#include "core/Configuration.h"
#include "../utils/math3d/primitives.h"
#include "math/MathUtility.h"
#include "math/MathUtility.h"
#include "applications/options/ExploreOptions.h"
#include "../RRTPlanner.h"
#include "../moves/CompositeMove.h"

using namespace std;

string TestDummy::name(){
	return "Rasmus' dummy test";
}

bool TestDummy::runTests(){
	torsionAccumTest();
    //sampleHIV1TAR();
    //sampleAndDisplay();
	//writeConfiguration();
	//printTorsions();
	//derivativeTest();

	//gsl_matrix* A = gsl_matrix_alloc(3,3);
	//gsl_matrix_set(A,0,0,0.5);gsl_matrix_set(A,0,1,1.0);gsl_matrix_set(A,0,2,0.0);
	//gsl_matrix_set(A,1,0,0.5);gsl_matrix_set(A,1,1,0.0);gsl_matrix_set(A,1,2,1.0);
	//gsl_matrix_set(A,2,0,0.5);gsl_matrix_set(A,2,1,0.0);gsl_matrix_set(A,2,2,0.0);
	//gsl_matrix_view B = gsl_matrix_submatrix(A, 0,1,3,2);
	//gsl_matrix* prod = gsl_matrix_alloc(3,3);

	//   gsl_blas_dgemm (CblasNoTrans, CblasTrans,
	//                     1.0, &B.matrix, &B.matrix,
	//                     0.0, prod);

	//gsl_matrix_cout(prod);

	//for(int i=0;i<6;i++){
	//	for(int j=0;j<10;j++){
	//		gsl_matrix_set(A,i,j,RandomN1P1());
	//	}
	//}
	//SVD svd(A);
	//cout<<"A:"<<endl;
	//gsl_matrix_cout(A);
	//svd.print();

	//gsl_matrix* A_dag = svd.pseudoInverse();
	//cout<<"A_dag:"<<endl;
	//gsl_matrix_cout(A_dag);
	//cout<<"A·A_dag:"<<endl;
	//gsl_matrix_cout(gsl_matrix_mul(A,A_dag));

	return true;
}

extern double COLLISION_FACTOR;

void TestDummy::torsionAccumTest(){
	SamplingOptions options;
	options.collisionFactor = 0.5;
	options.initialStructureFile = "/Users/rfonseca/init.pdb";
	options.hydrogenbondFile = "/Users/rfonseca/hbond.txt";
	options.extraCovBonds.push_back("123-628");
	options.maxRotation = 0.05;
	options.flexibleRibose = false;

    enableLogger("default");
    enableLogger("samplingStatus");
    COLLISION_FACTOR=options.collisionFactor;

    //string pdb_file = path + protein_name + ".pdb";
    string pdb_file = options.initialStructureFile;
	Protein* protein = new Protein();

	IO::readPdb( protein, pdb_file, options.extraCovBonds );
    log()<<"Molecule has "<<protein->Atom_list.size()<<" atoms\n";

	IO::readHbonds( protein, options.hydrogenbondFile );

	// Check for collision
	// This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
	protein->Initial_collisions = protein->getAllCollisions();

	IO::readRigidbody( protein );
	protein->buildRigidbodyTree(0, options.flexibleRibose);

    log("samplingStatus")<<"Molecule has:"<<endl;
    log("samplingStatus")<<"> "<<protein->Atom_list.size()<<" atoms"<<endl;
    log("samplingStatus")<<"> "<<protein->Initial_collisions.size()<<" initial collisions"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->DOF_num<<" DOFs of which "<<protein->m_spanning_tree->Cycle_DOF_num<<" are cycle-DOFs\n"<<endl;

    //Initiate planner
    //RRTPlanner planner( protein, options );
	cout<<"All good"<<endl;

	Configuration *pSmp = new Configuration(protein);
	protein->m_conf = pSmp;
	protein->m_conf_backup = pSmp;
	pSmp->m_parent = nullptr;

	string out_path = "./";
	string name = "init";

	//m_target = target;
	//Configuration *pTarget, *pClosestSmp, *pNewSmp;

	//Save initial file (for movie)
	string out_file = out_path + name + "_new_0.pdb";

	protein->SetConfiguration(pSmp);
	IO::writePdb(protein, out_file);

	Configuration* new_q = new Configuration(protein);
	new_q->m_ID = 0;
	new_q->m_parent = pSmp;
	new_q->m_tree_depth = pSmp->m_tree_depth+1;

	for(int it=1;it<2;it++){

		for (int i=0; i<new_q->m_DOF; ++i){
			new_q->m_f[i] = it*0.1745329252*4; //10deg
		}

		Configuration* pNewSmp = new_q;

		// Output sampled PDBs
		string out_file = out_path + name + "_new_" +
			Util::to_string(it) 
			+ ".pdb";

		IO::writePdb(pNewSmp->updatedProtein(), out_file);

		cout<< "> New structure: "<<name<<"_new_"<<it<<".pdb .. RMSD to initial: "<< setprecision(6)<<pNewSmp->m_rmsd_initial<<" .. Null-space dimension: "<<pNewSmp->CycleNullSpace->nullspaceSize<<endl<<endl;

	}



	exit(0);
}

void TestDummy::derivativeTest(){
	string path = "/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/Experiments/HIV1/free_1ANR/";
	string protein_name = "1ANR_3";
	string pdb_file = path + protein_name + ".pdb";
	vector<string> extraCovBonds;
	Protein* protein = new Protein();
	IO::readPdb( protein, pdb_file, extraCovBonds );

	string hbond_file = path + protein_name + ".pdb.out";
	IO::readHbonds_rnaview(protein, hbond_file, true);
	protein->Initial_collisions = protein->getAllCollisions();
	IO::readRigidbody( protein );
	protein->buildRigidbodyTree();
	//protein->BackupAtomPos();

	Configuration* conf = new Configuration(protein);

	//int dof = 20;//INteresting. Check residues
	for(int dof=0; dof<40; dof++)
	{
		bool wroteOrigin = false;
		double delta = 0.1;
		conf->m_f[dof] = delta;
		vector<Vector3> pose;//First-order approximation of delta change

		for(vector<Atom*>::iterator ait = protein->Atom_list.begin(); ait!=protein->Atom_list.end(); ++ait){
			Atom* atom = *ait;

			//Locate the RBVertex containing to atom
			RigidbodyGraphVertex* vertex = nullptr;
			map<unsigned int, RigidbodyGraphVertex*>::iterator vit;
			for(vit=protein->m_spanning_tree->Vertex_map.begin(); vit!=protein->m_spanning_tree->Vertex_map.end(); ++vit){
				RigidbodyGraphVertex* v = vit->second;
				if(v->Rb_ptr->containsAtom(atom)){
					vertex = v;
					break;
				}
			}


			//Trace from vertex to dof and put derivative in pose2
			while ( vertex != protein->m_spanning_tree->Root ) {
				RigidbodyGraphVertex* parent = vertex->Parent;
				Edge* p_edge = parent->Edges.find(vertex->VertexId)->second;

				if(parent->isRibose) {
					SugarVertex* v = reinterpret_cast<SugarVertex*>(parent);
					int dof_id = v->DOF_id;
					if(dof_id==dof){
						Vector3 derivative = v->computeJacobianEntry(p_edge, conf->m_f, atom->Position);
						derivative = atom->Position + (delta*derivative);
						pose.push_back( derivative );
						if(!wroteOrigin){
							cout<<"DOF "<<dof<<" corresponds to sugar residue "<<v->Rb_ptr->getAtom("C1'")->getResidue()->Id<<endl;
							wroteOrigin=true;
						}
						break;
					}
				}

				int dof_id = p_edge->DOF_id;
				if (dof_id==dof) {
					Atom* ea1 = p_edge->Bond->Atom1;
					Atom* ea2 = p_edge->Bond->Atom2;
					Vector3 derivative = MathUtility::ComputeJacobianEntry(ea1->Position,ea2->Position,atom->Position);
					derivative = atom->Position + (delta*derivative);
					pose.push_back( derivative );
					if(!wroteOrigin){
						cout<<"DOF "<<dof<<" corresponds to edge "<<p_edge->Bond<<endl;
						wroteOrigin = true;
					}
					break;
				}
				vertex = parent;
			}
			if(vertex==protein->m_spanning_tree->Root){
				pose.push_back(atom->Position);
			}

		}

		//		   cout<<"Edges.size : "<<protein->m_spanning_tree->Edges.size()<<endl;
		//		   for(vector<Edge*>::iterator eit = protein->m_spanning_tree->Edges.begin(); eit != protein->m_spanning_tree->Edges.end(); ++eit){
		//		   Edge* edge = *eit;
		////if(edge->DOF_id==dof){
		////	cout<<edge<<" has DOF_id "<<dof<<endl;
		////}
		////if(edge->Cycle_DOF_id==dof){
		////	cout<<edge<<" has Cycle_DOF_id "<<dof<<endl;
		////}
		//cout<<edge<<" has Cycle_DOF_id "<<edge->Cycle_DOF_id<<endl;
		//
		//}
		//map<unsigned int, RigidbodyGraphVertex*>::iterator vit;
		//for(vit=protein->m_spanning_tree->Vertex_map.begin(); vit!=protein->m_spanning_tree->Vertex_map.end(); ++vit){
		//RigidbodyGraphVertex* vert = vit->second;
		//if(vert->isRibose) {
		//SugarVertex* v = reinterpret_cast<SugarVertex*>(vert);
		////if(v->DOF_id==dof) cout<<"Sugarvertex "<<v->Rb_ptr->getAtom("C1'")->getResidue()->Id<<" has DOF_id "<<v->DOF_id<<endl;
		////if(v->Cycle_DOF_id==dof) cout<<"Sugarvertex "<<v->Rb_ptr->getAtom("C1'")->getResidue()->Id<<" has Cycle_DOF_id "<<v->Cycle_DOF_id<<endl;
		//cout<<"Sugarvertex "<<v->Rb_ptr->getAtom("C1'")->getResidue()->Id<<" has Cycle_DOF_id "<<v->Cycle_DOF_id<<endl;
		//}
		//}
		//


		protein->SetConfiguration(conf);

		char output_file [256];
		sprintf(output_file,"pose_dof_%d_real.pdb",dof);
		IO::writePdb(protein, output_file);

		sprintf(output_file,"pose_dof_%d_approx.pdb",dof);
		ofstream output(output_file);
		//for (vector<Atom*>::iterator atom_itr=protein->Atom_list.begin(); atom_itr!=protein->Atom_list.end(); ++atom_itr) {
		for(unsigned int a=0; a<protein->Atom_list.size(); a++){
			Atom* atom = protein->Atom_list[a];
			Vector3 pos = pose[a];
			Residue* res = atom->getResidue();

			char buffer[100];
			sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  ",
			        atom->getId(),atom->Original_name.c_str(),
			        res->Name.c_str(),res->getChainId().c_str(),res->getId(),
			        pos.x, pos.y, pos.z ,atom->getType().c_str());
			string line(buffer);
			output << line << endl;
		}
		output.close();
		conf->m_f[dof]=0.0;
		protein->SetConfiguration(conf);
	}

}

void TestDummy::printTorsions(){
	Vector3 p1( -2.10274, -8.64888, 5.80287);
	Vector3 p2( -0.596499, -8.49341, 5.26309);
	Vector3 p3( -0.745207, -8.68261, 3.68787);
	Vector3 p4( -0.16, -9.801, 3.024);

	cout<<TorsionalAngle(p1,p2,p3,p4)<<endl;

	p1.x = -2.103; p1.y = -8.649; p1.z = 5.803;
	p2.x = -0.596; p2.y = -8.493; p2.z = 5.263;
	p3.x = -0.745; p3.y = -8.683; p3.z = 3.688;
	cout<<TorsionalAngle(p1,p2,p3,p4)<<endl;
}

void printInterestingPart(Protein* p){
	Atom* C4_a = p->getAtom(6, "C4'"); log("test")<<C4_a<<endl;
	Atom* C3_a = p->getAtom(6, "C3'"); log("test")<<C3_a<<endl;
	Atom* C2_a = p->getAtom(6, "C2'"); log("test")<<C2_a<<endl;
	Atom* O3_a = p->getAtom(6, "O3'"); log("test")<<O3_a<<endl;
	Atom* P_a  = p->getAtom(7, "P");   log("test")<<P_a<<endl;
	Atom* O5_a = p->getAtom(7, "O5'"); log("test")<<O5_a<<endl;
	Atom* C5_a = p->getAtom(7, "C5'"); log("test")<<C5_a<<endl;
	Vector3 C4 = C4_a->Position;
	Vector3 C3 = C3_a->Position;
	Vector3 C2 = C2_a->Position;
	Vector3 O3 = O3_a->Position;
	Vector3 P  = P_a->Position;
	Vector3 O5 = O5_a->Position;
	Vector3 C5 = C5_a->Position;
	log("test")<<"C4-C3      : "<<(C4-C3).length()<<endl;
	log("test")<<"C3-O3      : "<<(C3-O3).length()<<endl;
	log("test")<<"O3-P       : "<<(O3-P).length()<<endl;
	log("test")<<"P-O5       : "<<(P-O5).length()<<endl;
	log("test")<<"O5-C5      : "<<(O5-C5).length()<<endl;
	log("test")<<"C4-C3-C2   : "<<VectorAngle(C4-C3, C2-C3)<<endl;
	log("test")<<"C4-C3-O3   : "<<VectorAngle(C4-C3, O3-C3)<<endl;
	log("test")<<"C3-O3-P    : "<<VectorAngle(C3-O3, P-O3)<<endl;
	log("test")<<"O3-P-O5    : "<<VectorAngle(O3-P, O5-P)<<endl;
	log("test")<<"P-O5-C5    : "<<VectorAngle(P-O5, C5-O5)<<endl;
	log("test")<<"C4-C3-O3-P : "<<TorsionalAngle(C4,C3,O3,P)<<endl;
	log("test")<<"C3-O3-P-O5 : "<<TorsionalAngle(C3,O3,P,O5)<<endl;
	log("test")<<"O3-P-O5-C5 : "<<TorsionalAngle(O3,P,O5,C5)<<endl;
}

void TestDummy::writeConfiguration(){
	string path = "/Users/rfonseca/Documents/RNAKGS/KGSX/branches/release-0.31-frozenRNA/tests/";
	string protein_name = "smallRNALoop";
	string pdb_file = path + protein_name + ".pdb";
	vector<string> extraCovBonds;
	Protein* protein = new Protein();
	IO::readPdb( protein, pdb_file, extraCovBonds );

	string hbond_output_file = path + "hbonds.in";
	IO::readHbonds(protein, hbond_output_file);
	protein->Initial_collisions = protein->getAllCollisions();
	IO::readRigidbody( protein );
	protein->buildRigidbodyTree();
	//protein->BackupAtomPos();

	//printInterestingPart(protein);

	Configuration* conf = new Configuration(protein);
	conf->m_f[2 ] = 0.4;
	//conf->m_f[5 ] = 0;//-2.85185 ;
	//conf->m_f[6 ] = 1 ;
	//conf->m_f[82] =  -0.0199798 ;
	//conf->m_f[83] =  -0.00661102 ;
	//conf->m_f[84] =  -0.0665515 ;
	//conf->m_f[85] =  -0.290697 ;
	//conf->m_f[86] =  0.16861 ;
	//conf->m_f[89] =  0.421758 ;
	//conf->m_f[90] =  0.295914 ;
	//conf->m_f[91] =  -0.283396 ;
	//conf->m_f[92] =  -0.080965 ;
	//conf->m_f[93] =  0;//0.183779 ; //C3'-O3'
	//conf->m_f[96] =  0;//-1.5708;//-0.110637 ; //O3'-P
	//conf->m_f[97] =  0;	//P-O5'
	enableLogger("debugRebuild");

	protein->SetConfiguration(conf);
	disableLogger("debugRebuild");
	//printInterestingPart(protein);

	//protein->resampleSugars(5,7, conf);

	IO::writePdb(protein, "confTest.pdb");
	log("test")<<"TestDummy::runTests: Wrote file confTest.pdb"<<endl;


}

extern double COLLISION_FACTOR;

void TestDummy::sampleAndDisplay(){
    srand(200);
    SamplingOptions options;
    options.initialStructureFile = "test_run/smallRNALoop2.pdb";

    options.hydrogenbondFile = "test_run/smallRNALoop2.pdb.out";
    options.hydrogenbondMethod = "rnaview";
    options.collisionFactor = 0.4;
    options.maxRotation = 0.08;
    options.rebuild_fragment_length = 2;
    options.rebuild_frequency = 0.10;
    options.rebuildAggression = 2;


    COLLISION_FACTOR=options.collisionFactor;

    //string pdb_file = path + protein_name + ".pdb";
    string pdb_file = options.initialStructureFile;
	vector<string> extraCovBonds;
	Protein* protein = new Protein();
	IO::readPdb( protein, pdb_file, extraCovBonds );
    log()<<"Molecule has "<<protein->Atom_list.size()<<" atoms\n";


    if(!options.annotationFile.empty())
        IO::readAnnotations(protein, options.annotationFile);

	if(options.hydrogenbondMethod=="user")
		IO::readHbonds( protein, options.hydrogenbondFile );
	else if(options.hydrogenbondMethod=="rnaview")
		IO::readHbonds_rnaview( protein, options.hydrogenbondFile, options.annotationFile.empty() );
	else if(options.hydrogenbondMethod=="dssr")
		IO::readHbonds_dssr( protein, options.hydrogenbondFile );

	// Check for collision
	// This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
	protein->Initial_collisions = protein->getAllCollisions();

	IO::readRigidbody( protein );
	protein->buildRigidbodyTree();

    log("samplingStatus")<<"Molecule has:"<<endl;
    log("samplingStatus")<<"> "<<protein->Atom_list.size()<<" atoms"<<endl;
    log("samplingStatus")<<"> "<<protein->Initial_collisions.size()<<" initial collisions"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->DOF_num<<" DOFs of which "<<protein->m_spanning_tree->Cycle_DOF_num<<" are cycle-DOFs\n";

	CompositeMove move(options);

    RRTPlanner planner( protein, options, move );
    //planner.Explore( options, target );
	//planner.m_protein->BackupAtomPos();

    Configuration *smp = *(planner.m_samples.begin()), *prevSmp;
	for(int i=0;i<100;i++){
        cout<<"Generating sample "<<i<<endl;
        prevSmp = smp;
		//smp = planner.Extend(smp, options);
		smp = move.performMove(smp, nullptr);
		if (smp) {
			string out_file = "test_run/output/smp_" +
			        Util::to_string(i) 
			        + ".pdb";
			IO::writePdb(smp->updatedProtein(), out_file);
		}else{
            smp = prevSmp;
            i--;
        }
    }

	delete protein;

    system("../../trunk/Scripts/combinePDBs.py test_run/smallRNALoop2.pdb test_run/output/smp_?.pdb test_run/output/smp_??.pdb > test_run/path.pdb");
    system("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL test_run/path.pdb -d 'intra_fit path'");
}


double minDonorAcceptorDistance2(Protein* p, int r1, int r2){
	string r1_donors[3], r1_acceptors[3];
	Atom* r1_C1 = p->getAtom(r1, "C1'");
	if(r1_C1->getResidue()->Name=="A"){
		r1_acceptors[0] = "N7";
		r1_acceptors[1] = "N1";
		r1_acceptors[2] = "N3";
		r1_donors[0] = "H61";
		r1_donors[1] = "H62";
		r1_donors[2] = "";
	}
	if(r1_C1->getResidue()->Name=="G"){
		r1_acceptors[0] = "N7";
		r1_acceptors[1] = "O6";
		r1_acceptors[2] = "N3";
		r1_donors[0] = "H1";
		r1_donors[1] = "H21";
		r1_donors[2] = "H22";
	}
	if(r1_C1->getResidue()->Name=="C"){
		r1_acceptors[0] = "N3";
		r1_acceptors[1] = "O2";
		r1_acceptors[2] = "";
		r1_donors[0] = "H41";
		r1_donors[1] = "H42";
		r1_donors[2] = "";
	}
	if(r1_C1->getResidue()->Name=="U"){
		r1_acceptors[0] = "O4";
		r1_acceptors[1] = "O2";
		r1_acceptors[2] = "";
		r1_donors[0] = "H3";
		r1_donors[1] = "";
		r1_donors[2] = "";
	}

	string r2_donors[3], r2_acceptors[3];
	Atom* r2_C1 = p->getAtom(r2, "C1'");
	if(r2_C1->getResidue()->Name=="A"){
		r2_acceptors[0] = "N7";
		r2_acceptors[1] = "N1";
		r2_acceptors[2] = "N3";
		r2_donors[0] = "H61";
		r2_donors[1] = "H62";
		r2_donors[2] = "";
	}
	if(r2_C1->getResidue()->Name=="G"){
		r2_acceptors[0] = "N7";
		r2_acceptors[1] = "O6";
		r2_acceptors[2] = "N3";
		r2_donors[0] = "H1";
		r2_donors[1] = "H21";
		r2_donors[2] = "H22";
	}
	if(r2_C1->getResidue()->Name=="C"){
		r2_acceptors[0] = "N3";
		r2_acceptors[1] = "O2";
		r2_acceptors[2] = "";
		r2_donors[0] = "H41";
		r2_donors[1] = "H42";
		r2_donors[2] = "";
	}
	if(r2_C1->getResidue()->Name=="U"){
		r2_acceptors[0] = "O4";
		r2_acceptors[1] = "O2";
		r2_acceptors[2] = "";
		r2_donors[0] = "H3";
		r2_donors[1] = "";
		r2_donors[2] = "";
	}

	double ret = 1000.0;
	for(int i=0;i<3&&r1_donors[i]!="";i++){
		Atom* a1 = p->getAtom(r1, r1_donors[i]);
		for(int j=0;j<3&&r2_acceptors[j]!="";j++){
			Atom* a2 = p->getAtom(r2, r2_acceptors[j]);
			double d = fabs( a1->Position.distance(a2->Position) - 2.0 );
			if(d<ret) ret=d;
		}
	}
	for(int i=0;i<3&&r1_acceptors[i]!="";i++){
		Atom* a1 = p->getAtom(r1, r1_acceptors[i]);
		for(int j=0;j<3&&r2_donors[j]!="";j++){
			Atom* a2 = p->getAtom(r2, r2_donors[j]);
			double d = fabs( a1->Position.distance(a2->Position) - 2.0 );
			if(d<ret) ret=d;
		}
	}
	return ret;
}


void TestDummy::sampleHIV1TAR(){
    srand(200);
    SamplingOptions options;
    options.initialStructureFile = "/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/Experiments/HIV1/free_1ANR/1ANR_3.pdb";
    options.hydrogenbondFile = 	"/Users/rfonseca/Documents/RNAKGS/KGSX/trunk/Experiments/HIV1/free_1ANR/hbonds_full_rebuildAll.txt";
    options.hydrogenbondMethod = "rnaview";
    options.collisionFactor = 0.5;
    options.maxRotation = 0.05;
    options.rebuild_fragment_length = 0;
    options.rebuild_frequency = 0.0;
    options.rebuildAggression = 1;


    COLLISION_FACTOR=options.collisionFactor;

    //string pdb_file = path + protein_name + ".pdb";
    string pdb_file = options.initialStructureFile;
	vector<string> extraCovBonds;
	Protein* protein = new Protein();
	IO::readPdb( protein, pdb_file, extraCovBonds );
    log()<<"Molecule has "<<protein->Atom_list.size()<<" atoms\n";


    if(!options.annotationFile.empty())
        IO::readAnnotations(protein, options.annotationFile);

	if(options.hydrogenbondMethod=="user")
		IO::readHbonds( protein, options.hydrogenbondFile );
	else if(options.hydrogenbondMethod=="rnaview")
		IO::readHbonds_rnaview( protein, options.hydrogenbondFile, options.annotationFile.empty() );
	else if(options.hydrogenbondMethod=="dssr")
		IO::readHbonds_dssr( protein, options.hydrogenbondFile );

	// Check for collision
	// This step is NECESSARY because it defines the original colliding atoms, and these atoms won't be considered as in collision during the sampling.
	protein->Initial_collisions = protein->getAllCollisions();

	IO::readRigidbody( protein );
	protein->buildRigidbodyTree();

    log("samplingStatus")<<"Molecule has:"<<endl;
    log("samplingStatus")<<"> "<<protein->Atom_list.size()<<" atoms"<<endl;
    log("samplingStatus")<<"> "<<protein->Initial_collisions.size()<<" initial collisions"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
    log("samplingStatus")<<"> "<<protein->m_spanning_tree->DOF_num<<" DOFs of which "<<protein->m_spanning_tree->Cycle_DOF_num<<" are cycle-DOFs\n";

	CompositeMove move(options);

	RRTPlanner planner( protein, options, move );
	//planner.Explore( options, target );
	//planner.m_protein->BackupAtomPos();

	Configuration *smp = *(planner.m_samples.begin()), *prevSmp;
	for(int i=0;i<100;i++){
		cout<<"Generating sample "<<i<<endl;
		prevSmp = smp;
		//smp = planner.Extend(smp, options);
		smp = move.performMove(smp, nullptr);
		if (smp) {

			//planner.m_protein->SetConfiguration(prevSmp);
			//double d1 = minDonorAcceptorDistance2(planner.m_protein, 30,35);
			//double d2 = minDonorAcceptorDistance2(planner.m_protein, 31,34);
			//double prevDistSum = d1+d2;
			//planner.m_protein->SetConfiguration(smp);
			//d1 = minDonorAcceptorDistance2(planner.m_protein, 30,35);
			//d2 = minDonorAcceptorDistance2(planner.m_protein, 31,34);
			//double distSum = d1+d2;
			//double distDiff = distSum - prevDistSum ;

			//if( Random01()<exp(-distDiff/0.5) ){
   //             cout<<"> Accepted: DistSum="<<distSum<<endl;
   //         }else{
			//	delete smp;
			//	smp = prevSmp;
			//	planner.m_protein->SetConfiguration(smp);
			//	i--;
			//	continue;
			//}

			string out_file = "test_run/output/smp_" +
			        Util::to_string(i) 
			        + ".pdb";
			IO::writePdb(smp->updatedProtein(), out_file);
		}else{
			smp = prevSmp;
			i--;
		}
	}

	delete protein;

	system("../../trunk/Scripts/combinePDBs.py /Users/rfonseca/Documents/RNAKGS/KGSX/trunk/Experiments/HIV1/free_1ANR/1ANR_3.pdb test_run/output/smp_?.pdb test_run/output/smp_{1..99}.pdb test_run/output/smp_{98..1}.pdb > test_run/path.pdb");
	system("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL test_run/path.pdb -d 'intra_fit path'");
}
