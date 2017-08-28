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
#include <iostream>

#include "Atom.h"
#include "Util.h"
#include "Bond.h"
#include "HBond.h"
#include "Rigidbody.h"

using namespace std;

Atom::Atom (const bool& hetatm, const string& name, const int& id, const Coordinate& pos, Residue* residue):
    m_hetatm(hetatm),
    m_name(name),
    m_id(id),
    m_position(pos),
    m_referencePosition(pos),
    m_parentResidue(residue),
    m_rigidbody(nullptr),
    m_biggerRigidbody(nullptr),
    m_bFactor(0)
{
//  On_sidechain = true;

  // Assign Element
  string::size_type char_start = m_name.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
  if ( char_start == string::npos ) {
    cerr << "Error: Empty atom name of atom " << m_id << endl;
    exit(1);
  }

//  if ( m_name.find_first_of("H")!=string::npos ) { // H atom
//    Element = atomH;
//  }
//  else if ( m_name.find("SE")!=string::npos ) { // SE atom
  if ( m_name.find("SE")!=string::npos ) { // SE atom
    m_element = atomSE;
  }
  else {
    char first_char = m_name[char_start];
    switch (first_char) {
      case 'C':
        m_element = atomC;
        break;
      case 'N':
        m_element = atomN;
        break;
      case 'O':
        m_element = atomO;
        break;
      case 'S':
        m_element = atomS;
        break;
      case 'H':
        m_element = atomH;
        break;
      default:
        m_element = atomOther;
        break;
    }
  }
}

const string& Atom::getName () const {
  return m_name;
}

int Atom::getId () const {
  return m_id;
}

double Atom::getMass () const {
  switch (m_element) {
    case atomC: return MASS_C;
    case atomH: return MASS_H;
    case atomN: return MASS_N;
    case atomO: return MASS_O;
    case atomS: return MASS_S;
    case atomSE:return MASS_SE;
    default:    return MASS_UNKNOWN;
  }
}

double Atom::getRadius() const {
  switch (m_element) {
    case atomC: return VDW_RADIUS_C;
    case atomH: return VDW_RADIUS_H;
    case atomN: return VDW_RADIUS_N;
    case atomO: return VDW_RADIUS_O;
    case atomS: return VDW_RADIUS_S;
    case atomSE:return VDW_RADIUS_SE;
    default:    return VDW_RADIUS_UNKNOWN;
  }
}

double Atom::getEpsilon() const {
  switch (m_element) {
    case atomC: return VDW_EPSILON_C;
    case atomH: return VDW_EPSILON_H;
    case atomN: return VDW_EPSILON_N;
    case atomO: return VDW_EPSILON_O;
    case atomS: return VDW_EPSILON_S;
    case atomSE:return VDW_EPSILON_SE;
    default:    return VDW_EPSILON_UNKNOWN;
  }
}

const std::string Atom::getElement() const{
  switch(m_element) {
    case atomC: return "C"; break;
    case atomH: return "H"; break;
    case atomN: return "N"; break;
    case atomO: return "O"; break;
    case atomS: return "S"; break;
    case atomSE:return "SE";break;
    default: return "?";
  }
}


void Atom::printSummaryInfo() const {
  cout << "Atom " << m_id << " " << m_name << " ";
  switch (m_element) {
    case atomC: cout << "C"; break;
    case atomH: cout << "H"; break;
    case atomN: cout << "N"; break;
    case atomO: cout << "O"; break;
    case atomS: cout << "S"; break;
    case atomSE:cout << "SE";break;
    default: cout << "UNKNOWN_ATOM";
  }

  cout << " Radius(" << getRadius() << ") " << (isBackboneAtom()?"Backbone":"Sidechain");
  cout << " Rigid("<<m_rigidbody->id() << ")";
  cout << " Pos(" << m_position.tostring() << ")" << endl;
  cout << "  " << Cov_bond_list.size() << " covalent bonds:";
  for (vector<Atom*>::const_iterator it=Cov_neighbor_list.begin(); it!=Cov_neighbor_list.end(); ++it) {
    Atom* neighbor = *it;
    cout << " (Res" << neighbor->getResidue()->getId() << " " << neighbor->getName() << " " << neighbor->getId() << ")";
  }
  cout << endl;
  cout << "  " << Hbond_list.size() << " h-bonds:";
  for (vector<Atom*>::const_iterator it=Hbond_neighbor_list.begin(); it!=Hbond_neighbor_list.end(); ++it) {
    Atom* neighbor = *it;
    cout << " (Res" << neighbor->getResidue()->getId() << " " << neighbor->getName() << " " << neighbor->getId() << ")";
  }
  cout << endl;
  cout << "  " << Second_cov_neighbor_list.size() << " 2nd immediate covalent neighbor atom IDs:";
  for (vector<Atom*>::const_iterator it=Second_cov_neighbor_list.begin(); it!=Second_cov_neighbor_list.end(); ++it)
    cout << " " << (*it)->m_id;
  cout << endl;
  //cout << "occ: " << occ << " u_iso " << u_iso << endl;
}

Atom* Atom::getBondNeighbor (Bond * bond) const {
  if ( bond->m_atom1 == this ) {
    return bond->m_atom2;
  }
  return bond->m_atom1;
}

Bond* Atom::getBond(Atom *other) const {
  for (vector<Bond*>::const_iterator bit=Cov_bond_list.begin(); bit!=Cov_bond_list.end(); ++bit) {
    Bond* bond = *bit;
    if ( bond->m_atom1 == this){
      if( bond->m_atom2 == other){
        return bond;
      }
    }
    else{
      if( bond->m_atom1 == other){
        return bond;
      }
    }
  }
  return nullptr;
}

Residue* Atom::getResidue () const{
  return m_parentResidue;
}

void Atom::addCovBond (Bond * cbond) {
  Cov_bond_list.push_back(cbond);
  Cov_neighbor_list.push_back( getBondNeighbor(cbond) );
}

void Atom::addHbond (Hbond * hbond) {
  Hbond_list.push_back(hbond);
  Hbond_neighbor_list.push_back( getBondNeighbor(hbond) );
}

void Atom::removeHbond (Hbond * hbond) {
  Atom* neighbor = getBondNeighbor(hbond);
  for (vector<Atom*>::iterator ait=Hbond_neighbor_list.begin(); ait!=Hbond_neighbor_list.end(); ++ait) {
    if ( (*ait)->m_id == neighbor->m_id ) {
      Hbond_neighbor_list.erase(ait);
      break;
    }
  }
  for (vector<Hbond *>::iterator hbit=Hbond_list.begin(); hbit != Hbond_list.end(); ++hbit) {
    if ( (*hbit)->isSame(hbond) ) {
      Hbond_list.erase(hbit);
      break;
    }
  }
}

double Atom::distanceTo (Atom* other) const {
  return m_position.distanceTo(other->m_position);
}

bool Atom::isWithinDistanceFrom (Atom* center, double dist) const {
  return m_position.isWithinSphere(center->m_position,dist);
}

//void Atom::setAsMainchainAtom () {
//  On_sidechain = false;
//}
//
//void Atom::setAsSidechainAtom () {
//  On_sidechain = true;
//}

Atom* Atom::getIthCovNeighbor (int i) const {
  if (Cov_neighbor_list.size()<i+1) {
    cerr << "Atom " << m_id <<"("<<this<<")"<< " has " << i << " covalent neighbors. "<<m_name << endl;
    exit(1);
  }
  return Cov_neighbor_list.at(i);
}

bool Atom::isHetatm() const {
  return m_hetatm;
}

bool Atom::isSidechainAtom () const {
//  return On_sidechain;
  return !isBackboneAtom();
}

void Atom::setRigidbody(Rigidbody* rb) {
  if(m_rigidbody!=nullptr)
    cerr<<"Atom::setRigidbody - OVERRIDING RIGID BODY NOW"<<endl;
  m_rigidbody = rb;
}

Rigidbody* Atom::getRigidbody() const{
  return m_rigidbody;
}

//void Atom::setBiggerRigidbody (Rigidbody* rb) {
//  m_biggerRigidbody = rb;
//}

//Rigidbody* Atom::getBiggerRigidbody() const {
//  return m_biggerRigidbody;
//}

string Atom::getType () const {
  return getName().substr(0,1);
}

bool Atom::isCovNeighbor (Atom* other) const {
  for (vector<Atom*>::const_iterator it=Cov_neighbor_list.begin(); it!=Cov_neighbor_list.end(); ++it) {
    if ( (*it) == other )
      return true;
  }
  return false;
}

bool Atom::isHbondNeighbor (Atom* other) const {
  for (vector<Atom*>::const_iterator it=Hbond_neighbor_list.begin(); it!=Hbond_neighbor_list.end(); ++it) {
    if ( (*it) == other )
      return true;
  }
  return false;
}

bool Atom::isSecondCovNeighbor (Atom* other) const {
  for (vector<Atom*>::const_iterator it=Second_cov_neighbor_list.begin(); it!=Second_cov_neighbor_list.end(); ++it) {
    if ( (*it) == other )
      return true;
  }
  return false;
}

Atom* Atom::getFirstCovNeighbor () const {
  return getIthCovNeighbor(0);
}

bool Atom::compare (Atom* atom1, Atom* atom2) {
  return (atom1->m_id < atom2->m_id);
}

bool Atom::compareName (string name) const {
  if( getName() == name )
    return true;
  return false;
}

bool Atom::compareType (string type) const {
  if( getType() == type )
    return true;
  return false;
}

//vector<Rigidbody*> Atom::getSameRigidbody (Atom* another) const {
//  vector<Rigidbody*> rbs;
//  for (vector<Rigidbody*>::const_iterator it1=Rigidbody_list.begin(); it1!=Rigidbody_list.end(); ++it1) {
//    for (vector<Rigidbody*>::const_iterator it2=another->Rigidbody_list.begin(); it2!=another->Rigidbody_list.end(); ++it2) {
//      if ( (*it1) == (*it2) )
//        rbs.push_back(*it1);
//    }
//  }
//  return rbs;
//}

void Atom::setBFactor(float val) {
  m_bFactor = val;
}

float Atom::getBFactor() {
  return m_bFactor;
}

void Atom::setOccupancy(float val) {
  m_occupancy = val;
}

float Atom::getOccupancy() {
  return m_occupancy;
}

bool Atom::inSameRigidbody (Atom* another) const {
  if(getRigidbody()==another->getRigidbody()) return true;
  for(auto const& bond: getRigidbody()->m_bonds){
    if(bond->m_atom1==another) return true;
    if(bond->m_atom2==another) return true;
  }
  return false;

  //for (vector<Rigidbody*>::const_iterator it1=Rigidbody_list.begin(); it1!=Rigidbody_list.end(); ++it1) {
  //  for (vector<Rigidbody*>::const_iterator it2=another->Rigidbody_list.begin(); it2!=another->Rigidbody_list.end(); ++it2) {
  //    if ( (*it1) == (*it2) )
  //      return true;
  //  }
  //}
  //return false;
}


bool Atom::isBackboneAtom () const {
  return ( getName() == "C" || getName() == "CA" || getName() == "N" || getName() == "O" );
}

bool Atom::isHeavyAtom () const {
  return (getType() != "H");
}


bool Atom::isCollisionCheckAtom (string collisionCheckAtoms ) const {
  if( collisionCheckAtoms == "backbone" )
    return isBackboneAtom();
  else if( collisionCheckAtoms == "heavy" )
    return isHeavyAtom();
  else if( collisionCheckAtoms == "all" )
    return true;
  else if( collisionCheckAtoms == "none" )
    return false;
  else
    return true;
}

/* Return heavy-atom valence for this atom */
int Atom::getHAV() const {
  int ret = 0;
  for (auto const &neighbor: Cov_neighbor_list) {
    if (neighbor->isHeavyAtom())
      ++ret;
  }
  return ret;
}

std::vector<Atom*> Atom::heavyAtomNeighbors() const{
  std::vector<Atom*> heavyNeighbors;
  for(auto const& neighbor: Cov_neighbor_list){
    if (neighbor->isHeavyAtom())
      heavyNeighbors.push_back(neighbor);
  }
  return heavyNeighbors;
}

/* Output operator overloading for easy printing. */

ostream& operator<<(ostream& os, const Atom& a) {
  //os<<"Atom["<<a.m_id<<","<<a.m_parentResidue->m_name<<"_"<<a.m_parentResidue->m_id<<"_"<<a.m_name<<"]";
  os<<a.getResidue()->getName()<<"_"<<a.getResidue()->getId()<<"_"<<a.getName();
  return os;
}

ostream& operator<<(ostream& os, const Atom* a) {
  return os<<*a;
}

