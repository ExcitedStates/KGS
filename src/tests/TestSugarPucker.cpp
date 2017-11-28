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

#include "TestSugarPucker.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>
#include <sstream>


#include "../IO.h"
#include "core/Molecule.h"
#include "core/HBond.h"
#include "core/Atom.h"
#include "../Logger.h"

using namespace std;

/** Verify that RNA ribose and protein proline works ok. This implies that:
  *  1)  C4' is entry => Positions 0-5 = (C4', C3', C2', C1', O4')
  *  2)  C3' is entry => Positions 0-5 = (C3', C4', O4', C1', C2')
  *  3)  C1' is entry => Positions 0-5 = (C1', C2', C3', C4', O4')
  *  4)  CA is entry => Positions 0-5 = (CA, N, CG, CD, CB)
  *  5)  N is entry => Positions 0-5 = (N, CA, CB, CD, CG)
  *  6)  Only position 2 and 3 and attached moves.
  *  7)  tau = 0  => |tau_3-Am| ≈ 0
  *  8)  tau = pi => |tau_3+Am| ≈ 0
  *  9)  (tau = pi/2 v tau = 3pi/2) => |tau_3| ≈ 0
  * 10)  C1'(pi/2) ≠ C1'(3pi/2)
  * 11)  C1'(pi/2) ≈ C1'(5pi/2)  (Periodicity)
  * TODO: Test position derivatives.
  * TODO: Test transformation propagation
  * TODO: Test what happens if ribose/proline is root
  */
bool TestSugarPucker::runTests(){
	bool passed = true;
    if(testRNA()) log("test")<<left<<setw(60)<<"TestSugarPucker::testRNA:"<<"passed"<<endl;
    else { log("test")<<left<<setw(60)<<"TestSugarPucker::testRNA:"<<"failed"<<endl; passed = false;}
    if(testProtein()) log("test")<<left<<setw(60)<<"TestSugarPucker::testProtein:"<<"passed"<<endl;
    else { log("test")<<left<<setw(60)<<"TestSugarPucker::testProtein:"<<"failed"<<endl; passed = false;}
    if(testTransformationPropagation()) 	log("test")<<left<<setw(60)<<"TestSugarPucker::testTransformationPropagation:"<<"passed"<<endl;
    else { log("test")<<left<<setw(60)<<"TestSugarPucker::testTransformationPropagation:"<<"failed"<<endl; passed = false;}
    return passed;
}
bool TestSugarPucker::testProtein(){
    string pdb_file = "tests/polypro.pdb";
	vector<string> extraCovBonds;
    Protein* protein = new Protein();
    IO::readPdb( protein, pdb_file, extraCovBonds );


    IO::readRigidbody( protein );
    protein->buildRigidbodyTree(14);
    //protein->m_spanning_tree->printForSpringy();


    SugarVertex* v1 = reinterpret_cast<SugarVertex*>(protein->m_spanning_tree->getVertex(12)); // CA is entry
    SugarVertex* v2 = reinterpret_cast<SugarVertex*>(protein->m_spanning_tree->getVertex(16));  // N is entry

    Atom* a_1_0 = v1->ringAttached[0].front();
    Atom* a_1_1 = v1->ringAttached[1].front();
    Atom* a_1_2 = v1->ringAttached[2].front();
    Atom* a_1_3 = v1->ringAttached[3].front();
    Atom* a_1_4 = v1->ringAttached[4].front();

    Atom* a_2_0 = v2->ringAttached[0].front();
    Atom* a_2_1 = v2->ringAttached[1].front();
    Atom* a_2_2 = v2->ringAttached[2].front();
    Atom* a_2_3 = v2->ringAttached[3].front();
    Atom* a_2_4 = v2->ringAttached[4].front();

    // 4)  CA is entry => Positions 0-5 = (CA, N, CD, CG, CB)
    if(a_1_0->Name!="CA"){ log("test")<<"TestSugarPucker::testProtein: Test 4) expected CA but was "<<a_1_0->Name<<endl; return false; }
    if(a_1_1->Name!="N" ){ log("test")<<"TestSugarPucker::testProtein: Test 4) expected N but was " <<a_1_1->Name<<endl; return false; }
    if(a_1_2->Name!="CD"){ log("test")<<"TestSugarPucker::testProtein: Test 4) expected CD but was "<<a_1_2->Name<<endl; return false; }
    if(a_1_3->Name!="CG"){ log("test")<<"TestSugarPucker::testProtein: Test 4) expected CG but was "<<a_1_3->Name<<endl; return false; }
    if(a_1_4->Name!="CB"){ log("test")<<"TestSugarPucker::testProtein: Test 4) expected CB but was "<<a_1_4->Name<<endl; return false; }

    //5)  N is entry => Positions 0-5 = (N, CA, CB, CD, CG)
    if(a_2_0->Name!="N" ){ log("test")<<"TestSugarPucker::testProtein: Test 5) expected N but was "<<a_2_0->Name<<endl; return false; }
    if(a_2_1->Name!="CA"){ log("test")<<"TestSugarPucker::testProtein: Test 5) expected CA but was "<<a_2_1->Name<<endl; return false; }
    if(a_2_2->Name!="CB"){ log("test")<<"TestSugarPucker::testProtein: Test 5) expected CB but was "<<a_2_2->Name<<endl; return false; }
    if(a_2_3->Name!="CG"){ log("test")<<"TestSugarPucker::testProtein: Test 5) expected CG but was "<<a_2_3->Name<<endl; return false; }
    if(a_2_4->Name!="CD"){ log("test")<<"TestSugarPucker::testProtein: Test 5) expected CD but was "<<a_2_4->Name<<endl; return false; }


//    CConfiguration* conf = new CConfiguration(protein->m_spanning_tree->DOF_num);
    Configuration* conf = new Configuration(protein);//----------------CHANGED---------------
    Vector3 p_1_0 = a_1_0->Position;
    Vector3 p_1_1 = a_1_1->Position;
    Vector3 p_1_2 = a_1_2->Position;
    Vector3 p_1_3 = a_1_3->Position;
    Vector3 p_1_4 = a_1_4->Position;

    Vector3 p_2_0 = a_2_0->Position;
    Vector3 p_2_1 = a_2_1->Position;
    Vector3 p_2_2 = a_2_2->Position;
    Vector3 p_2_3 = a_2_3->Position;
    Vector3 p_2_4 = a_2_4->Position;

    conf->m_f[1]+=0.31415;
    protein->SetConfiguration(conf);

    // 6)  Only position 2 and 3 moves.
    if(p_1_0.distance(a_1_0->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_0<<" to be "<<a_1_0->Position<<endl; return false; }
    if(p_1_1.distance(a_1_1->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_1<<" to be "<<a_1_1->Position<<endl; return false; }
    if(p_1_2.distance(a_1_2->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_2<<" to change"<<endl; return false; }
    if(p_1_3.distance(a_1_3->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_3<<" to change"<<endl; return false; }
    if(p_1_4.distance(a_1_4->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_4<<" to be "<<a_1_4->Position<<endl; return false; }

    conf->m_f[1]=0;
    conf->m_f[3]+=0.31415;
    protein->SetConfiguration(conf);
    if(p_2_0.distance(a_2_0->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_0<<" to be "<<a_2_0->Position<<endl; return false; }
    if(p_2_1.distance(a_2_1->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_1<<" to be "<<a_2_1->Position<<endl; return false; }
    if(p_2_2.distance(a_2_2->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_2<<" to change"<<endl; return false; }
    if(p_2_3.distance(a_2_3->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_3<<" to change"<<endl; return false; }
    if(p_2_4.distance(a_2_4->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_4<<" to be "<<a_2_4->Position<<endl; return false; }

    // 7)  tau = 0  => |tau_3-Am| ≈ 0
    double initTau = acos(v2->initTorsion/v2->Amplitude)*v2->initUp;
    conf->m_f[3]=-initTau;
    protein->SetConfiguration(conf);
    double tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3-v2->Amplitude)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be "<<v2->Amplitude<<" but was "<<tau_3<<endl;return false; }

    // 8)  tau = pi  => |tau_3+Am| ≈ 0
    conf->m_f[3]=-initTau+M_PI;
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3+v2->Amplitude)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be "<<v2->Amplitude<<" but was "<<tau_3<<endl;return false; }

    // 9)  (tau = pi/2 v tau = 3pi/2) => |tau_3| ≈ 0
    conf->m_f[3]=-initTau+(M_PI/2);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp1 = a_2_3->Position;

    conf->m_f[3]=-initTau+(3*M_PI/2.0);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp2 = a_2_3->Position;

    //10)  C1'(pi/2) ≠ C1'(3pi/2)
    if(C1_tmp1.distance(C1_tmp2)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 10) expected C1' to move"<<endl;return false; }

    //11)  C1'(pi/2) ≈ C1'(5pi/2)  (Periodicity)
    conf->m_f[3]=-initTau+(5*M_PI/2.0);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp3 = a_2_3->Position;
    if(C1_tmp1.distance(C1_tmp3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 11) expected C1' to be periodic"<<endl;return false; }

    return true;
}

bool TestSugarPucker::testRNA(){
    string pdb_file = "tests/smallRNALoop.pdb";
	vector<string> extraCovBonds;
    Protein* protein = new Protein();
    IO::readPdb( protein, pdb_file, extraCovBonds );

    //Manually add a h-bond
    Atom* hatom = protein->getAtom(25);
    Atom* donor = hatom->getFirstCovNeighbor();
    Atom* oatom = protein->getAtom(174);
    Atom* AA = oatom->getFirstCovNeighbor();
    ProteinHbond* new_hb = new ProteinHbond(hatom,oatom,donor,AA);
    protein->addHbond(new_hb);

    IO::readRigidbody( protein );
    protein->buildRigidbodyTree(18);

    SugarVertex* v1 = reinterpret_cast<SugarVertex*>(protein->m_spanning_tree->getVertex(22)); // C4' is entry
    SugarVertex* v2 = reinterpret_cast<SugarVertex*>(protein->m_spanning_tree->getVertex(8));  // C3' is entry
    SugarVertex* v3 = reinterpret_cast<SugarVertex*>(protein->m_spanning_tree->getVertex(15)); // C1' is entry

    Atom* a_1_0 = v1->ringAttached[0].front();
    Atom* a_1_1 = v1->ringAttached[1].front();
    Atom* a_1_2 = v1->ringAttached[2].front();
    Atom* a_1_3 = v1->ringAttached[3].front();
    Atom* a_1_4 = v1->ringAttached[4].front();

    Atom* a_2_0 = v2->ringAttached[0].front();
    Atom* a_2_1 = v2->ringAttached[1].front();
    Atom* a_2_2 = v2->ringAttached[2].front();
    Atom* a_2_3 = v2->ringAttached[3].front();
    Atom* a_2_4 = v2->ringAttached[4].front();

    Atom* a_3_0 = v3->ringAttached[0].front();
    Atom* a_3_1 = v3->ringAttached[1].front();
    Atom* a_3_2 = v3->ringAttached[2].front();
    Atom* a_3_3 = v3->ringAttached[3].front();
    Atom* a_3_4 = v3->ringAttached[4].front();

    // 1)  C4' is entry => Positions 0-5 = (C4', C3', C2', C1', O4')
    if(a_1_0->Name!="C4'"){ log("test")<<"TestSugarPucker::testRNA: Test 1) expected C4' but was "<<a_1_0->Name<<endl; return false; }
    if(a_1_1->Name!="C3'"){ log("test")<<"TestSugarPucker::testRNA: Test 1) expected C3' but was "<<a_1_1->Name<<endl; return false; }
    if(a_1_2->Name!="C2'"){ log("test")<<"TestSugarPucker::testRNA: Test 1) expected C2' but was "<<a_1_2->Name<<endl; return false; }
    if(a_1_3->Name!="C1'"){ log("test")<<"TestSugarPucker::testRNA: Test 1) expected C1' but was "<<a_1_3->Name<<endl; return false; }
    if(a_1_4->Name!="O4'"){ log("test")<<"TestSugarPucker::testRNA: Test 1) expected O4' but was "<<a_1_4->Name<<endl; return false; }

    // 2)  C3' is entry => Positions 0-5 = (C3', C4', O4', C1', C2')
    if(a_2_0->Name!="C3'"){ log("test")<<"TestSugarPucker::testRNA: Test 2) expected C3' but was "<<a_2_0->Name<<endl; return false; }
    if(a_2_1->Name!="C4'"){ log("test")<<"TestSugarPucker::testRNA: Test 2) expected C4' but was "<<a_2_1->Name<<endl; return false; }
    if(a_2_2->Name!="O4'"){ log("test")<<"TestSugarPucker::testRNA: Test 2) expected O4' but was "<<a_2_2->Name<<endl; return false; }
    if(a_2_3->Name!="C1'"){ log("test")<<"TestSugarPucker::testRNA: Test 2) expected C1' but was "<<a_2_3->Name<<endl; return false; }
    if(a_2_4->Name!="C2'"){ log("test")<<"TestSugarPucker::testRNA: Test 2) expected C2' but was "<<a_2_4->Name<<endl; return false; }

    // 3)  C1' is entry => Positions 0-5 = (C1', C2', C3', C4', O4')
    if(a_3_0->Name!="C1'"){ log("test")<<"TestSugarPucker::testRNA: Test 3) expected C1' but was "<<a_3_0->Name<<endl; return false; }
    if(a_3_1->Name!="C2'"){ log("test")<<"TestSugarPucker::testRNA: Test 3) expected C2' but was "<<a_3_1->Name<<endl; return false; }
    if(a_3_2->Name!="C3'"){ log("test")<<"TestSugarPucker::testRNA: Test 3) expected C3' but was "<<a_3_2->Name<<endl; return false; }
    if(a_3_3->Name!="C4'"){ log("test")<<"TestSugarPucker::testRNA: Test 3) expected C4' but was "<<a_3_3->Name<<endl; return false; }
    if(a_3_4->Name!="O4'"){ log("test")<<"TestSugarPucker::testRNA: Test 3) expected O4' but was "<<a_3_4->Name<<endl; return false; }

//    CConfiguration* conf = new CConfiguration(protein->m_spanning_tree->DOF_num);
    Configuration* conf = new Configuration(protein);//----------------CHANGED-------------
    Vector3 p_1_0 = a_1_0->Position;
    Vector3 p_1_1 = a_1_1->Position;
    Vector3 p_1_2 = a_1_2->Position;
    Vector3 p_1_3 = a_1_3->Position;
    Vector3 p_1_4 = a_1_4->Position;

    Vector3 p_2_0 = a_2_0->Position;
    Vector3 p_2_1 = a_2_1->Position;
    Vector3 p_2_2 = a_2_2->Position;
    Vector3 p_2_3 = a_2_3->Position;
    Vector3 p_2_4 = a_2_4->Position;

    conf->m_f[3]+=0.31415;
    protein->SetConfiguration(conf);

    // 6)  Only position 2 and 3 moves.
    if(p_1_0.distance(a_1_0->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_0<<" to be "<<a_1_0->Position<<endl; return false; }
    if(p_1_1.distance(a_1_1->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_1<<" to be "<<a_1_1->Position<<endl; return false; }
    if(p_1_2.distance(a_1_2->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_2<<" to change"<<endl; return false; }
    if(p_1_3.distance(a_1_3->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_3<<" to change"<<endl; return false; }
    if(p_1_4.distance(a_1_4->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_1_4<<" to be "<<a_1_4->Position<<endl; return false; }

    conf->m_f[3]=0;
    conf->m_f[1]+=0.31415;
    protein->SetConfiguration(conf);
    if(p_2_0.distance(a_2_0->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_0<<" to be "<<a_2_0->Position<<endl; return false; }
    if(p_2_1.distance(a_2_1->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_1<<" to be "<<a_2_1->Position<<endl; return false; }
    if(p_2_2.distance(a_2_2->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_2<<" to change"<<endl; return false; }
    if(p_2_3.distance(a_2_3->Position)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_3<<" to change"<<endl; return false; }
    if(p_2_4.distance(a_2_4->Position)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 6) expected "<<p_2_4<<" to be "<<a_2_4->Position<<endl; return false; }

    // 7)  tau = 0  => |tau_3-Am| ≈ 0
    double initTau = acos(v2->initTorsion/v2->Amplitude)*v2->initUp;
    conf->m_f[1]=-initTau;
    protein->SetConfiguration(conf);
    double tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3-v2->Amplitude)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be "<<v2->Amplitude<<" but was "<<tau_3<<endl;return false; }

    // 8)  tau = pi  => |tau_3+Am| ≈ 0
    conf->m_f[1]=-initTau+M_PI;
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3+v2->Amplitude)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be "<<v2->Amplitude<<" but was "<<tau_3<<endl;return false; }

    // 9)  (tau = pi/2 v tau = 3pi/2) => |tau_3| ≈ 0
    conf->m_f[1]=-initTau+(M_PI/2);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp1 = a_2_3->Position;

    conf->m_f[1]=-initTau+(3*M_PI/2.0);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp2 = a_2_3->Position;

    //10)  C1'(pi/2) ≠ C1'(3pi/2)
    if(C1_tmp1.distance(C1_tmp2)<0.001){ log("test")<<"TestSugarPucker::testRNA: Test 10) expected C1' to move"<<endl;return false; }

    //11)  C1'(pi/2) ≈ C1'(5pi/2)  (Periodicity)
    conf->m_f[1]=-initTau+(5*M_PI/2.0);
    protein->SetConfiguration(conf);
    tau_3 = TorsionalAngle(a_2_4->Position, a_2_0->Position, a_2_1->Position, a_2_2->Position);
    if(fabs(tau_3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Expected tau_3 to be 0 but was "<<tau_3<<endl;return false; }
    Vector3 C1_tmp3 = a_2_3->Position;
    if(C1_tmp1.distance(C1_tmp3)>0.001){ log("test")<<"TestSugarPucker::testRNA: Test 11) expected C1' to be periodic"<<endl;return false; }

    return true;
}

void fillMap(map<string,double>& map, Protein* prot){
	//Fill in bond-lengths
	for(vector<Edge*>::iterator eit = prot->m_spanning_tree->Edges.begin(); eit != prot->m_spanning_tree->Edges.end(); eit++){
		Edge* e = *eit;

		stringstream lengthKey;
		lengthKey<<e->Bond->Atom1->getName()<<"_"<<e->Bond->Atom1->getId()<<"-"<<e->Bond->Atom1->getName()<<"_"<<e->Bond->Atom1->getId();
		map[lengthKey.str()] = (e->Bond->Atom1->Position-e->Bond->Atom2->Position).length();

		stringstream torsionKey;
		Atom* a1 = e->Bond->Atom1->Cov_neighbor_list[0]==e->Bond->Atom2?e->Bond->Atom1->Cov_neighbor_list[1]:e->Bond->Atom1->Cov_neighbor_list[0];
		Atom* a2 = e->Bond->Atom1;
		Atom* a3 = e->Bond->Atom2;
		Atom* a4 = e->Bond->Atom2->Cov_neighbor_list[0]==e->Bond->Atom1?e->Bond->Atom2->Cov_neighbor_list[1]:e->Bond->Atom2->Cov_neighbor_list[0];

		//Make sure a1 is the covalent neighbor of Atom1 with lexicographically smallest name
		for(vector<Atom*>::iterator ait = e->Bond->Atom1->Cov_neighbor_list.begin(); ait!=e->Bond->Atom1->Cov_neighbor_list.end(); ait++){
			Atom* a = *ait;
			if(a!=e->Bond->Atom2 && a->Name<a1->Name) a1 = a;
		}
		//Make sure a4 is the covalent neighbor of Atom2 with lexicographically smallest name
		for(vector<Atom*>::iterator ait = e->Bond->Atom2->Cov_neighbor_list.begin(); ait!=e->Bond->Atom2->Cov_neighbor_list.end(); ait++){
			Atom* a = *ait;
			if(a!=e->Bond->Atom1 && a->Name<a4->Name) a4 = a;
		}
		torsionKey<<a1->Name<<"_"<<a1->getId()<<"-";
		torsionKey<<a2->Name<<"_"<<a2->getId()<<"-";
		torsionKey<<a3->Name<<"_"<<a3->getId()<<"-";
		torsionKey<<a4->Name<<"_"<<a4->getId();
		//torsionKey<<e->DOF_id;
		map[torsionKey.str()] = e->Bond->getTorsion();


		stringstream angleKey1;
		angleKey1<<a1->Name<<"_"<<a1->getId()<<"-";
		angleKey1<<a2->Name<<"_"<<a2->getId()<<"-";
		angleKey1<<a3->Name<<"_"<<a3->getId();
		map[angleKey1.str()] = VectorAngle(a1->Position-a2->Position, a3->Position-a2->Position);

		stringstream angleKey2;
		angleKey1<<a2->Name<<"_"<<a2->getId()<<"-";
		angleKey1<<a3->Name<<"_"<<a3->getId()<<"-";
		angleKey1<<a4->Name<<"_"<<a4->getId();
		map[angleKey2.str()] = VectorAngle(a2->Position-a3->Position, a4->Position-a3->Position);
	}
}

void findMismatches(map<string,double>& map1, map<string,double>& map2, vector<string>& ret){
	for(map<string,double>::iterator mit = map1.begin(); mit!=map1.end(); mit++){
		string key1 = mit->first;
		double value1 = mit->second;
		double value2 = map2[key1];
		if(fabs(value1-value2)>0.001){
			ret.push_back(key1);
		}
	}
}

bool TestSugarPucker::testTransformationPropagation(){
    string pdb_file = "tests/smallRNALoop.pdb";
	vector<string> extraCovBonds;
    Protein* protein = new Protein();
    IO::readPdb( protein, pdb_file, extraCovBonds );

    //Manually add a h-bond
    Atom* hatom = protein->getAtom(25);
    Atom* donor = hatom->getFirstCovNeighbor();
    Atom* oatom = protein->getAtom(174);
    Atom* AA = oatom->getFirstCovNeighbor();
    ProteinHbond* new_hb = new ProteinHbond(hatom,oatom,donor,AA);
    protein->addHbond(new_hb);

    IO::readRigidbody( protein );
    protein->buildRigidbodyTree(0);

	//protein->m_spanning_tree->printForSpringy();
	IO::writePdb(protein, "confTest_0.pdb");

	map<string, double> initialMap;
	map<string, double> finalMap;
	vector<string> mismatches;
	fillMap(initialMap, protein);

//    CConfiguration* conf = new CConfiguration(protein->m_spanning_tree->DOF_num);
    Configuration* conf = new Configuration(protein);//---------------CHANGED-------------
	conf->m_f[3] = 1;//Ribose conformation in residue 4
	for(int i=10; i<conf->m_DOF; i++)
		conf->m_f[i] = 1;//Ribose conformation in residue 4
	//conf->m_f[28]+=1;log("test")<<"changing 31: ";
	for(vector<Edge*>::iterator eit = protein->m_spanning_tree->Edges.begin(); eit != protein->m_spanning_tree->Edges.end(); eit++){
		Edge* e = *eit;
		if(e->DOF_id==28) log("test")<<e->Bond<<endl;
	}
	//conf->m_f[32]+=1;
	protein->SetConfiguration(conf);
	IO::writePdb(protein, "confTest.pdb");
	fillMap(finalMap, protein);
	findMismatches(initialMap, finalMap, mismatches);
	if(!mismatches.empty()){
		bool passed = true;
		for(vector<string>::iterator mit = mismatches.begin(); mit!=mismatches.end(); mit++){
			string key = *mit;
			double newVal = initialMap[key]+1; if(newVal>CTK_PI) newVal-=2*CTK_PI;
			if(fabs( newVal-finalMap[key] )>0.0001){
				log("test")<<"TransformationPropagation 1) "<<key<<" should be "<<initialMap[key]<<" but is "<<finalMap[key]<<endl;
				passed = false;
			}
		}
		return passed;
	}

	return true;
}
