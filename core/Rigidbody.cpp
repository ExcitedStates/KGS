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
#include "Rigidbody.h"

using namespace std;

#include <iostream>

//---------------------------------------------------------
// Constructors and Destructors
Rigidbody::Rigidbody () {m_rbVertex = nullptr;}
Rigidbody::Rigidbody (unsigned int id) : rbId_(id) {m_rbVertex = nullptr;}
Rigidbody::~Rigidbody () {}
//---------------------------------------------------------
// Accessor id
unsigned int Rigidbody::id () const {
	return rbId_;
}
//---------------------------------------------------------
// Mutator id
void Rigidbody::id (unsigned int id) {
	rbId_ = id;
}

void Rigidbody::addAtom (Atom* atom) {
	Atoms.push_back(atom);
	atom->setRigidbody(this);
}

void Rigidbody::addBond (Bond * bond) {
	Bonds.push_back(bond);
}

void Rigidbody::setVertex (KinVertex* vertex) {
  m_rbVertex = vertex;
}

KinVertex* Rigidbody::getVertex (){
  return m_rbVertex;
}

void Rigidbody::makeBiggerRigidBody (Rigidbody* rb, Bond * bond) {
	//combine two rigid bodies that are linked by a constrained bond
	//bond is the constrained linking bond between the two bodies
	for( vector<Atom*>::iterator ait = rb->Atoms.begin(); ait != rb->Atoms.end(); ait++){
		if( ! this->containsAtom((*ait)) ){
			Atoms.push_back((*ait));
			(*ait)->setBiggerRigidbody(this);
		}
		//else{//already in our body, but also in the other one, so remove it!
		//	//	This happens for atoms of connecting bonds!
		//	(*ait)->removeBiggerRigidbody(rb);
		//}
	}
	for(vector<Bond *>::iterator bit = rb->Bonds.begin(); bit != rb->Bonds.end(); bit++){
		if( (*bit)->Atom1->getId() != bond->Atom1->getId() && !(*bit)->isHbond() ){
			Bonds.push_back((*bit));
		}
	}
}


void Rigidbody::makeBiggerRigidBody (Rigidbody* rb) {
	// forms a bigger rigid body, starting with this one as a base
	// therefore, no checking necessary
	for( vector<Atom*>::iterator ait = rb->Atoms.begin(); ait != rb->Atoms.end(); ait++){
		if( ! this->containsAtom((*ait)) ){
			Atoms.push_back((*ait));
			(*ait)->setBiggerRigidbody(this);
		}
		//else{//already in our body, but also in the other one, so remove it!
//	//		cout<<"Already in here in body "<<this->id()<<endl;
		//	(*ait)->removeBiggerRigidbody(rb);
		//}
	}

	for(vector<Bond *>::iterator bit = rb->Bonds.begin(); bit != rb->Bonds.end(); bit++){
			Bonds.push_back((*bit));
	}
}

int Rigidbody::size () const {
	return Atoms.size();
}

void Rigidbody::print () const {
	cout << "RigidBody ID " << id() << ": ";
	//cout << "RigidBody " << RbId << ": ";
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
			cout << (*it)->getId() << " ";
	}
	cout << endl;
}

void Rigidbody::printAtoms () const {
	cout << "\t RigidBody ID " << id() << ": " << endl;
	cout << "\t\t :" << (isMainchainRb()?"isOnMainchain":"isNOTOnMainchain") << endl;
	//cout << "\t\t :" << (isMainchainRb_?"isOnMainchain":"isNOTOnMainchain") << endl;
	cout << "\t\t numAtoms: " << size() << endl;
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it)
		cout << "\t\t Atom: " << (*it)->getResidue()->getId() << " " << (*it)->getId() << " " << (*it)->getName() << endl;
}

void Rigidbody::printAtomsBonds () const {
	cout << "\t RigidBody ID " << id() << ": " << endl;
	cout << "\t\t :" << (isMainchainRb()?"isOnMainchain":"isNOTOnMainchain") << endl;
	cout << "\t\t numAtoms: " << size() << endl;
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it)
		cout << "\t\t Atom: " << (*it)->getResidue()->getId() << " " << (*it)->getId() << " " << (*it)->getName() << endl;

	cout << "\t\t numBonds = " << Bonds.size() << endl;
	for (vector<Bond *>::const_iterator it=Bonds.begin(); it != Bonds.end(); ++it) {
        	cout << "\t\t Bond: " << (*it)->Atom1->getResidue()->getId() << " " << (*it)->Atom1->getId() << " " << (*it)->Atom1->getName() <<
					    "\t ---------> \t" 
				   << (*it)->Atom2->getResidue()->getId() << " " << (*it)->Atom2->getId() << " " << (*it)->Atom2->getName();
		cout << "\t\t Type: " << ((*it)->BondType=="HB"?"HB":"CV") << endl;
	}
}

bool Rigidbody::containsResidue( Residue* res ) const {
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
		if( (*it)->getResidue()==res )
			return true;
	}
	return false;
}

bool Rigidbody::containsAtom (Atom* atom) const {
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
		if ( (*it)==atom )
			return true;
	}
	return false;
}

//bool Rigidbody::containsAtomAtPosition( const clipper::Coord_orth& pos ) const {
//	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it)
//		if ( pos[0] >= (*it)->m_Position[0] - (*it)->Vdw_radius / 100 &&
//		     pos[0] <= (*it)->m_Position[0] + (*it)->Vdw_radius / 100 &&
//		     pos[1] >= (*it)->m_Position[1] - (*it)->Vdw_radius / 100 &&
//		     pos[1] <= (*it)->m_Position[1] + (*it)->Vdw_radius / 100 &&
//		     pos[2] >= (*it)->m_Position[2] - (*it)->Vdw_radius / 100 &&
//		     pos[2] <= (*it)->m_Position[2] + (*it)->Vdw_radius / 100 )
//			return true;
//	return false;
//}

bool Rigidbody::containsMainchainAtoms() const {
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
		if( (*it)->getName() == "C" || (*it)->getName() == "CA" || (*it)->getName() == "N" || (*it)->getName() == "O" )
			return true;
	}
	return false;
}

bool Rigidbody::containsAlongMainchainAtoms() const {
	for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
		if( (*it)->getName() == "C" || (*it)->getName() == "CA" || (*it)->getName() == "N" )
			return true;
	}
	return false;
}

bool Rigidbody::isWithinResidueRange( unsigned int resid1, unsigned int resid2 ) const {
	if ( resid2>=resid1 )
		for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
			if( (*it)->getResidue()->getId() >= resid1 && (*it)->getResidue()->getId() <= resid2 )
				return true;
		}
	return false;
}

bool Rigidbody::isWithinTwoResidueRanges( unsigned int resid1, unsigned int resid2,
					  unsigned int resid3, unsigned int resid4 ) const {
	if ( resid2>=resid1 && resid4>=resid3 )
		for (vector<Atom*>::const_iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
			if( ( (*it)->getResidue()->getId() >= resid1 && (*it)->getResidue()->getId() <= resid2 )
			 || ( (*it)->getResidue()->getId() >= resid3 && (*it)->getResidue()->getId() <= resid4 ) )
				return true;
		}
	return false;
}

// Assume at least one of C and N AND one CA must be in the Rigidbody
void Rigidbody::setMainchainRb () {
	int count1 = 0;
	int count2 = 0;
	for (vector<Atom*>::iterator it=Atoms.begin(); it!=Atoms.end(); ++it) {
		if( (*it)->getName() == "C" || (*it)->getName() == "N" )
			count1++ ;
		if( (*it)->getName() == "CA" )
			count2++ ;
	}
	if ( count1 > 0 && count2 > 0 )
		isMainchainRb_ = true;
	else
		isMainchainRb_ = false;
}

bool Rigidbody::isMainchainRb () const {
	return isMainchainRb_;
}


Atom* Rigidbody::getAtom(string name){
	for(int i=0;i<Atoms.size();i++){
		if(Atoms[i]->getName()==name)
			return Atoms[i];
	}
	return nullptr;
}


ostream& operator<<(ostream& os, const Rigidbody& rb){
	os<<"Rigidbody[";
    //int count = 0;
    for(vector<Atom*>::const_iterator it = rb.Atoms.begin(); it!=rb.Atoms.end(); ++it){
        os<<*(*it)<<",";
        //if(++count>2) { os<<"..,"; break; }
    }
	os<<"ID: "<<rb.id()<<"]";
	return os;
}

ostream& operator<<(ostream& os, const Rigidbody* rb){
	os<<*rb;
	return os;
}
