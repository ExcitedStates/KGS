
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <sstream>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "core/Atom.h"
#include "JacobianRelated.h"
#include "Util.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "RunFirst.h"
#include "math/MathUtility.h"
#include "DisjointSets.h"
#include "core/ProteinHBond.h"
#include "Logger.h"
#include "Selection.h"


using namespace std;

double sugar_up(Vector3& p0, Vector3& p1, Vector3& p2, Vector3& p3, Vector3& p4 ){

	Vector3 n = (p2-p1)/(p2-p1).length();
	Vector3 S1 = p2+n*( cos(CTK_PI-VectorAngle(p3-p2, p1-p2))*(p2-p3).length() );
	double S2off = n.dot(S1-p4);
	Vector3 S2 = p4+( n*S2off );

	double p3side = cross( p3-S1,S2-S1 ).dot(n);
	return Signum(p3side);
}

bool sugar_valid(Vector3& p0, Vector3& p1, Vector3& p2, Vector3& p3, Vector3& p4, double Am){
	double torsion = Am;

	double angles[5];
	double dists[5];
	angles[0] = VectorAngle(p1-p0,p4-p0);
	angles[1] = VectorAngle(p2-p1,p0-p1);
	angles[2] = VectorAngle(p3-p2,p1-p2);
	angles[3] = VectorAngle(p4-p3,p2-p3);
	angles[4] = VectorAngle(p0-p4,p3-p4);
	dists[0] = (p0-p1).length();
	dists[1] = (p1-p2).length();
	dists[2] = (p2-p3).length();
	dists[3] = (p3-p4).length();
	dists[4] = (p4-p0).length();

	Vector3 p0x = p1-p0; normalize(p0x);
	Vector3 p0z = cross(p0x, p4-p0); normalize(p0z);
	Vector3 p0y = cross(p0z, p0x); 
	Vector3 p2_ = p1+(cos(CTK_PI-angles[1])*dists[1])*p0x + (sin(CTK_PI-angles[1])*cos(torsion)*dists[1])*p0y + (sin(CTK_PI-angles[1])*sin(torsion)*dists[1])*p0z;

	Vector3 n = (p2_-p1)/dists[1];
	Vector3 S1 = p2_+n*( cos(CTK_PI-angles[2])*dists[2] );
	double S2off = n.dot(S1-p4);
	Vector3 S2 = p4 + n*S2off;
	double l = (S1-S2).length();
	double r1 = sin(CTK_PI-angles[2])*dists[2];
	double r2 = ( dists[3]*dists[3] - S2off*S2off );//Still need to take the square root
	if( r2<0 ) return false;
	r2 = sqrt(r2);
	if(r1+r2<l) return false;
	return true;
}

double sugar_amplitude(Vector3& p0, Vector3& p1, Vector3& p2, Vector3& p3, Vector3& p4){
	double lower = 0, upper = 3.14;
	while(upper-lower>0.0001){
		double tmpAm = (lower+upper)/2;
		if(sugar_valid(p0,p1,p2,p3,p4,tmpAm)) 	{lower = tmpAm;}
		else									{upper = tmpAm;}
	}
	return lower;
}

int main( int argc, char* argv[] ){
    enableLogger("sugar");

    log("sugar")<<"#    m_protein\t";
    log("sugar")<<"     residue\t";
    log("sugar")<<"       delta\t";
    log("sugar")<<"        tau0\t";
    log("sugar")<<"        tau1\t";
    log("sugar")<<"        tau2\t";
    log("sugar")<<"        tau3\t";
    log("sugar")<<"        tau4\t";
    log("sugar")<<"         chi\t";
    log("sugar")<<" O4'-C1'-C2'\t";
    log("sugar")<<" C1'-C2'-C3'\t";
    log("sugar")<<" C2'-C3'-C4'\t";
    log("sugar")<<" C3'-C4'-O4'\t";
    log("sugar")<<" C4'-O4'-C1'\t";
    log("sugar")<<"     C1'-C2'\t";
    log("sugar")<<"     C2'-C3'\t";
    log("sugar")<<"     C3'-C4'\t";
    log("sugar")<<"     C4'-O4'\t";
    log("sugar")<<"     O4'-C1'\t";
    log("sugar")<<"        tau3\t";
    log("sugar")<<"          up\t";
    log("sugar")<<"          Am\t";
    log("sugar")<<"         tau\t";
    log("sugar")<<endl;

	int start = 1;	
	vector<string> extraCovBonds;
	if(argc>3 && string(argv[1])=="--extraCovBonds") {
		Selection::split( argv[2],",", extraCovBonds );
		start = 3;
	}

	for(int a=start;a<argc;a++){
		char* tmp = realpath(argv[a], nullptr);
		if(tmp==nullptr){ cerr<<argv[a]<<" is not a valid PDB-file"<<endl; exit(-1); }
		string pdb_file(tmp);
		int nameSplit = pdb_file.find_last_of("/\\");
		string path = pdb_file.substr(0,nameSplit+1);
		string protein_name = pdb_file.substr(nameSplit+1);
		int pos = protein_name.rfind(".pdb");
		if(pos!=string::npos) protein_name = protein_name.substr(0,pos);

		Molecule * protein = new Molecule();
        IO::readPdb( protein, argv[a], extraCovBonds );

		IO::readRigidbody( protein );
		protein->buildSpanningTree();
		
		for(auto vit=protein->m_spanning_tree->Vertex_map.begin(); vit!=protein->m_spanning_tree->Vertex_map.end(); ++vit){
			KinVertex* v = vit->second;
			if(v->isRibose){
				SugarVertex* vs = reinterpret_cast<SugarVertex*>(v);
                //log("sugar")<<"Sugr"<<setw(14)<<left<<vs->DOF_id<<"\t";
				Vector3 C1 = v->m_rigidbody->getAtom("C1'")->m_Position;
				Vector3 C2 = v->m_rigidbody->getAtom("C2'")->m_Position;
				Vector3 C3 = v->m_rigidbody->getAtom("C3'")->m_Position;
				Vector3 C4 = v->m_rigidbody->getAtom("C4'")->m_Position;
				Vector3 O4 = v->m_rigidbody->getAtom("O4'")->m_Position;
				Vector3 C5 = v->m_rigidbody->getAtom("C5'")->m_Position;
				Vector3 O3 = v->m_rigidbody->getAtom("O3'")->m_Position;
				Atom* scNAtom = v->m_rigidbody->getAtom("N9");
				Atom* scCAtom = nullptr;
				if(scNAtom!=nullptr){//purine
					for(int i=1;i<scNAtom->Cov_neighbor_list.size();i++) 
						if(scNAtom->Cov_neighbor_list[i]->getName()=="C4"){
							scCAtom = scNAtom->Cov_neighbor_list[i];
							break;
						}
				}else{//pyrimidine
					scNAtom = v->m_rigidbody->getAtom("N1");
					for(int i=1;i<scNAtom->Cov_neighbor_list.size();i++) 
						if(scNAtom->Cov_neighbor_list[i]->getName()=="C2"){
							scCAtom = scNAtom->Cov_neighbor_list[i];
							break;
						}
				}
				Vector3 scN = scNAtom->m_Position;
				Vector3 scC = scCAtom->m_Position;
				//double up = sugar_up(C4,C3,C2,C1,O4);
				//double tau3 = TorsionalAngle(O4, C4, C3, C2);
				//double Am = sugar_amplitude(C4,C3,C2,C1,O4);
				//double tau = acos(tau3/Am)*up;

                log("sugar")<<setw(12)<<protein_name<<"\t"; //1
                log("sugar")<<setw(12)<<v->m_rigidbody->Atoms[0]->getResidue()->getId()<<"\t"; //2
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(C5, C4, C3, O3)<<"\t";//delta .. 3
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(C4, O4, C1, C2)<<"\t";//tau0 .. 4
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(O4, C1, C2, C3)<<"\t";//tau1 .. 5
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(C1, C2, C3, C4)<<"\t";//tau2 .. 6
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(C2, C3, C4, O4)<<"\t";//tau3 .. 7
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(C3, C4, O4, C1)<<"\t";//tau4 .. 8
                log("sugar")<<setw(12)<<setprecision(4)<<TorsionalAngle(O4, C1, scN, scC)<<"\t";//chi .. 9
                log("sugar")<<setw(12)<<setprecision(4)<<VectorAngle(O4-C1, C2-C1)<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<VectorAngle(C1-C2, C3-C2)<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<VectorAngle(C2-C3, C4-C3)<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<VectorAngle(C3-C4, O4-C4)<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<VectorAngle(C4-O4, C1-O4)<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<(C1-C2).length()<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<(C2-C3).length()<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<(C3-C4).length()<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<(C4-O4).length()<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<(O4-C1).length()<<"\t";
                log("sugar")<<setw(12)<<setprecision(4)<<vs->initTorsion<<"\t"; 			//tau3 .. 20
                log("sugar")<<setw(12)<<setprecision(4)<<vs->initUp<<"\t";					//UP   .. 21
                log("sugar")<<setw(12)<<setprecision(4)<<vs->Amplitude<<"\t";				//Am   .. 22
                log("sugar")<<setw(12)<<setprecision(4)<<acos(vs->initTorsion/vs->Amplitude)*vs->initUp<<"\t";
                //log("sugar")<<setw(12)<<setprecision(4)<<tau3<<"\t";//vs->initTorsion<<"\t"; //tau3 .. 20
                //log("sugar")<<setw(12)<<setprecision(4)<<up<<"\t";//vs->initUp<<"\t";
                //log("sugar")<<setw(12)<<setprecision(4)<<Am<<"\t";//vs->Amplitude<<"\t";
                //log("sugar")<<setw(12)<<setprecision(4)<<tau;
                log("sugar")<<endl;


			}
		}

	}
}
