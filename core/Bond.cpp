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
#include "Bond.h"
#include "Atom.h"
#include "math/MathUtility.h"
#include "Logger.h"

#include <iostream>

using namespace std;

Bond::Bond(Atom* atom1, Atom* atom2, std::string bond_type) {
	if (atom1->getId()<=atom2->getId()) {
		Atom1 = atom1;
		Atom2 = atom2;
	}
	else {
		Atom1 = atom2;
		Atom2 = atom1;
	}
	BondType = bond_type;
	Bars = 5; // default. Cannot compute the exact number of bars when the bond is just created because it doesn't know the number of covalent neighbors yet.
	constrained = false; ///only true if the bond is locked due to constraints (like hbonds)
}

Bond::Bond(Bond & bond) {
	Atom1 = bond.Atom1;
	Atom2 = bond.Atom2;
	BondType = bond.BondType;
	Bars = bond.Bars;
	constrained = bond.constrained;
}

Bond::Bond() {
	Atom1 = 0;
	Atom2 = 0;
	BondType = "UNDEFINED";
	Bars = 0;
	constrained=false;
}

Bond::~Bond() {}

void Bond::print () {
	cout << "Bond(" << Atom1->getId() << "," << Atom2->getId() << "," << BondType << "," << Bars << ")";
}

bool Bond::isLocked () {
	bool result = false;
	if ( Atom1->m_element==atomC && Atom1->Cov_neighbor_list.size()<=3 && Atom2->m_element==atomC && Atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	else if ( Atom1->m_element==atomC && Atom1->Cov_neighbor_list.size()<=3 && Atom2->m_element==atomN && Atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	else if ( Atom1->m_element==atomN && Atom1->Cov_neighbor_list.size()<=3 && Atom2->m_element==atomC && Atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	if(constrained)
		result = true;

	return result;
}

bool Bond::isPeptideBond () {
	if ( (Atom1->getName()=="C" && Atom2->getName()=="N") || (Atom1->getName()=="N" && Atom2->getName()=="C") )
		return true;
	return false;
}

bool Bond::isHbond () {
	return BondType=="HB";
}

double Bond::getTorsion() {
    
	int atom_id1 = Atom1->getId(); // due to the assertion of Bond, atom_id1 must be smaller than atom_id2
	int atom_id2 = Atom2->getId();
	Atom* atom3 = nullptr;
	for (vector<Atom*>::iterator aitr=Atom1->Cov_neighbor_list.begin(); aitr!=Atom1->Cov_neighbor_list.end(); ++aitr) {
		if ( (*aitr)->getId()==atom_id2 ) continue;
		if ( atom3==nullptr || (*aitr)->getId()<atom3->getId() ) {
			atom3 = *aitr;
		}
	}
	Atom* atom4 = nullptr;
	for (vector<Atom*>::iterator aitr=Atom2->Cov_neighbor_list.begin(); aitr!=Atom2->Cov_neighbor_list.end(); ++aitr) {
		if ( (*aitr)->getId()==atom_id1 ) continue;
		if ( atom4==nullptr || (*aitr)->getId()<atom4->getId() ) {
			atom4 = *aitr;
		}
	}

	double ret = 0.0;
	if(atom3 != nullptr && atom4 != nullptr){//only measure it if four covalently bonded atoms exist
		ret = TorsionalAngle(atom3->m_position,Atom1->m_position,Atom2->m_position,atom4->m_position); // in radians
		ret = formatRangeRadian(ret);
	}
	else{
//		log("dominik")<<"Setting global torsion to zero at bond between atom "<<Atom1->getId()<<" and "<<Atom2->getId()<<endl;
	}
	return ret;
    
}

ostream& operator<<(ostream& os, const Bond & b) {
	return os<<"Bond["<<b.Atom1<<"-"<<b.Atom2<<"]";
}

ostream& operator<<(ostream& os, const Bond * b) {
	return os<<*b;
}
