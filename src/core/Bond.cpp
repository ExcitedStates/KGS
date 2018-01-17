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


#include "Bond.h"
#include "Atom.h"
#include "math/MathUtility.h"
#include "Logger.h"

#include <iostream>

using namespace std;

Bond::Bond(Atom* atom1, Atom* atom2, std::string bond_type) {
	if (atom1->getId()<=atom2->getId()) {
		m_atom1 = atom1;
		m_atom2 = atom2;
	}
	else {
		m_atom1 = atom2;
		m_atom2 = atom1;
	}
	m_bondType = bond_type;
	m_bars = 5; // default. Cannot compute the exact number of bars when the bond is just created because it doesn't know the number of covalent neighbors yet.
	rigidified = false; ///only true if the bond is locked due to constraints (like hbonds)
}

Bond::Bond(Bond & bond) {
	m_atom1 = bond.m_atom1;
	m_atom2 = bond.m_atom2;
	m_bondType = bond.m_bondType;
	m_bars = bond.m_bars;
	rigidified = bond.rigidified;
}

Bond::Bond() {
	m_atom1 = 0;
	m_atom2 = 0;
	m_bondType = "UNDEFINED";
	m_bars = 0;
	rigidified=false;
}

Bond::~Bond() {}

void Bond::print () {
//	cout << "Bond(" << Atom1->getId() << "," << m_atom2->getId() << "," << m_bondType << "," << m_bars << ")";
	cout << "Bond(" << m_atom1 << "," << m_atom2 << "," << m_bondType << "," << m_bars << ")";
}

bool Bond::isLocked() const{
	bool result = false;
	if ( m_atom1->m_element==atomC && m_atom1->Cov_neighbor_list.size()<=3 && m_atom2->m_element==atomC && m_atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	else if ( m_atom1->m_element==atomC && m_atom1->Cov_neighbor_list.size()<=3 && m_atom2->m_element==atomN && m_atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	else if ( m_atom1->m_element==atomN && m_atom1->Cov_neighbor_list.size()<=3 && m_atom2->m_element==atomC && m_atom2->Cov_neighbor_list.size()<=3 )
		result = true;
	if(rigidified)
		result = true;

	return result;
}

bool Bond::isPeptideBond() const {
	if ( (m_atom1->getName()=="C" && m_atom2->getName()=="N") || (m_atom1->getName()=="N" && m_atom2->getName()=="C") )
		return true;
	return false;
}

bool Bond::isHBond() const{
	return m_bondType=="HB";
}
bool Bond::isDBond() const{
	return m_bondType=="DB";
}
bool Bond::isHydrophobicBond() const{
	return m_bondType=="HYB";
}


double Bond::getTorsion() {
    
	int atom_id1 = m_atom1->getId(); // due to the assertion of Bond, atom_id1 must be smaller than atom_id2
	int atom_id2 = m_atom2->getId();
	Atom* atom3 = nullptr;
	for (vector<Atom*>::iterator aitr=m_atom1->Cov_neighbor_list.begin(); aitr!=m_atom1->Cov_neighbor_list.end(); ++aitr) {
		if ( (*aitr)->getId()==atom_id2 ) continue;
		if ( atom3==nullptr || (*aitr)->getId()<atom3->getId() ) {
			atom3 = *aitr;
		}
	}
	Atom* atom4 = nullptr;
	for (vector<Atom*>::iterator aitr=m_atom2->Cov_neighbor_list.begin(); aitr!=m_atom2->Cov_neighbor_list.end(); ++aitr) {
		if ( (*aitr)->getId()==atom_id1 ) continue;
		if ( atom4==nullptr || (*aitr)->getId()<atom4->getId() ) {
			atom4 = *aitr;
		}
	}

	double ret = 0.0;
	if(atom3 != nullptr && atom4 != nullptr){//only measure it if four covalently bonded atoms exist
		ret = TorsionalAngle(atom3->m_position,m_atom1->m_position,m_atom2->m_position,atom4->m_position); // in radians
		ret = formatRangeRadian(ret);
	}
	else{
//		log("planner")<<"Setting global torsion to zero at bond between atom "<<Atom1->getId()<<" and "<<m_atom2->getId()<<endl;
	}
	return ret;
    
}

///Compare IDs of two bonds, used to sort them, lowest ID goes first
bool Bond::compareIDs(Bond *bond1, Bond *bond2) {
	//Determine min/max IDs in case bonds were not from min to max ID (happens for hbonds)
	int minID1, maxID1, minID2, maxID2;
	minID1 = bond1->m_atom1->getId();
	if (bond1->m_atom2->getId() < minID1 ){
		maxID1 = minID1;
		minID1 = bond1->m_atom2->getId();
	}
	else
		maxID1 = bond1->m_atom2->getId();

	minID2 = bond2->m_atom1->getId();
	if (bond2->m_atom2->getId() < minID2 ){
		maxID2 = minID2;
		minID2 = bond2->m_atom2->getId();
	}
	else
		maxID2 = bond2->m_atom2->getId();
	//Sort
	if(minID1 < minID2 )
		return true;
	if(minID1 > minID2)
		return false;
	return maxID1 < maxID2;
}

ostream& operator<<(ostream& os, const Bond & b) {
	return os<<"Bond["<<b.m_atom1<<"-"<<b.m_atom2<<"]";
}

ostream& operator<<(ostream& os, const Bond * b) {
	return os<<*b;
}
