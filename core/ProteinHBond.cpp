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
#include "ProteinHBond.h"
#include "core/Atom.h"
#include "math/MathUtility.h"

using namespace Math3D;
using namespace std;

Hbond::Hbond(Atom* hatom, Atom* acceptor, Atom* donor, Atom* aa, double energy) {
	Hatom = hatom;
	Acceptor = acceptor;
	Donor = donor;
	AA = aa;
	Energy = energy;

	// super-class attributes
	Atom1 = hatom;
	Atom2 = acceptor;
	BondType = "HB";
	Bars = 5;
	constrained = false;
	iniLength = VectorLength(Atom1->m_Position,Atom2->m_Position);
	iniOrientLeft = getLeftAngle();
	iniOrientRight = getRightAngle();

    Vector3 x,y,z;
    coordinateSystem(hatom, x,y,z);
    Vector3 d = acceptor->m_Position -hatom->m_Position;
    idealA[0] = x.dot(d);
    idealA[1] = y.dot(d);
    idealA[2] = z.dot(d);
    coordinateSystem(acceptor, x,y,z);
    idealH[0] = -x.dot(d);
    idealH[1] = -y.dot(d);
    idealH[2] = -z.dot(d);
}

Hbond::Hbond(Hbond & hbond) {
	Hatom = hbond.Hatom;
	Acceptor = hbond.Acceptor;
	Donor = hbond.Donor;
	AA = hbond.AA;
	Energy = hbond.Energy;

	// super-class attributes
	Atom1 = hbond.Atom1;
	Atom2 = hbond.Atom2;
	BondType = hbond.BondType;
	Bars = hbond.Bars;
	constrained = hbond.constrained;
	iniLength = hbond.iniLength;
	iniOrientLeft = hbond.iniOrientLeft;
	iniOrientRight = hbond.iniOrientRight;

	idealA = hbond.idealA;
	idealH = hbond.idealH;
}

bool Hbond::isSame (Hbond * b2) {
	if ( Hatom->getName() == b2->Hatom->getName() &&
		Hatom->getResidue()->getId() == b2->Hatom->getResidue()->getId() &&
		Donor->getName()==b2->Donor->getName() &&
		Donor->getResidue()->getId() == b2->Donor->getResidue()->getId() &&
		Acceptor->getName()==b2->Acceptor->getName() &&
		Acceptor->getResidue()->getId() == b2->Acceptor->getResidue()->getId() &&
		AA->getName()==b2->AA->getResidue()->getName() &&
		AA->getResidue()->getId() == b2->AA->getResidue()->getId() ){
			return true;
	}
	return false;
}

Vector3 Hbond::getIdealHPoint(){
    Vector3 x,y,z;
    coordinateSystem(Acceptor, x,y,z);
    return Acceptor->m_Position + (x*idealH[0]) + (y*idealH[1]) + (z*idealH[2]);
}

Vector3 Hbond::getIdealAcceptorPoint(){
    Vector3 x,y,z;
    coordinateSystem(Hatom, x,y,z);
    return Hatom->m_Position + (x*idealA[0]) + (y*idealA[1]) + (z*idealA[2]);
}

void Hbond::coordinateSystem(Atom* a, Vector3& x, Vector3& y, Vector3& z ){
	Atom* a1 = a;
    Atom* a2 = a->Cov_neighbor_list[0];
    Atom* a3 = nullptr;
    if( a1->Cov_neighbor_list.size()>1) a3 = a1->Cov_neighbor_list[1];
    else{
        if(a2->Cov_neighbor_list[0]!=a1)
            a3 = a2->Cov_neighbor_list[0];
        else
            a3 = a2->Cov_neighbor_list[1];
    }

    x = a1->m_Position -a2->m_Position;
    y = a3->m_Position -a2->m_Position;
    z = cross(x,y);
}

double Hbond::getLeftAngle() {
    Atom *a1, *a2, *a3;
    a1 = Bond::Atom1->Cov_neighbor_list[0] == Atom2 ? Atom1->Cov_neighbor_list[1] : Atom1->Cov_neighbor_list[0];
    a2 = Bond::Atom1;
    a3 = Bond::Atom2;
    //Make sure a1 is the covalent neighbor of Atom1 with lexicographically smallest name
    for(vector<Atom*>::iterator ait = a2->Cov_neighbor_list.begin(); ait!=a2->Cov_neighbor_list.end(); ait++){
        Atom* a = *ait;
        if(a!=a3 && a->getName()<a1->getName()) a1 = a;
    }
    double ret=VectorAngle(a2->m_Position - a1->m_Position, a3->m_Position - a2->m_Position);
    return ret;
}

double Hbond::getRightAngle() {
    Atom *a1, *a2, *a3;
    a1 = Bond::Atom1;
    a2 = Bond::Atom2;
    a3 = Bond::Atom2->Cov_neighbor_list[0] == Atom1 ? Atom2->Cov_neighbor_list[1] : Atom2->Cov_neighbor_list[0];
    //Make sure a3 is the covalent neighbor of Atom2 with lexicographically smallest name
    for(vector<Atom*>::iterator ait = a2->Cov_neighbor_list.begin(); ait!=a2->Cov_neighbor_list.end(); ait++){
        Atom* a = *ait;
        if(a!=a1 && a->getName()<a3->getName()) a3 = a;
    }
    double ret=VectorAngle(a2->m_Position - a1->m_Position, a3->m_Position - a2->m_Position);
    return ret;
}

