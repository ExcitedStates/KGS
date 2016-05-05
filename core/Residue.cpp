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

#include "Residue.h"
#include "Util.h"
#include "core/Atom.h"
#include "core/Chain.h"

using namespace std;

Residue::Residue(const string& name, const int& id, const Chain* parent_chain, const string& sse_type):
    Name(name),
    Id(id),
    m_parentChain(parent_chain),
    SSE_type(sse_type)
{
  m_prevResidue = nullptr;
  m_nextResidue = nullptr;
}

Residue::~Residue() {
}

void Residue::setLastResidue(Residue* last) {
	m_prevResidue = last;
	if (last!=nullptr) {
		last->m_nextResidue = this;
	}
}

void Residue::setNextResidue(Residue* next) {
	m_nextResidue = next;
	if (next!=nullptr) {
		next->m_prevResidue = this;
	}
}

void Residue::printSummaryInfo() const {
	cout << "Residue " << Name << " " << Id << " " << SSE_type << " " << sidechainSize() << endl;
//	for (map<string,Atom*>::const_iterator it= name_to_atom_map.begin(); it != name_to_atom_map.end(); ++it) {
//  it->second->printSummaryInfo();
  for (auto const& atom: m_atoms) {
		atom->printSummaryInfo();
	}
}

void Residue::print() const {
	cout << "Residue " << Name << " " << Id << " " << SSE_type << " " << sidechainSize() << endl;
}


Atom* Residue::addAtom(const std::string& atomName,
                      const int& atomId,
                      const Coordinate& position)
{
  Atom* atom = new Atom(atomName, atomId, position, this);
	//name_to_atom_map[atom->getName()] = atom;
  m_atoms.push_back(atom);
  return atom;
}

int Residue::getId () const {
	return Id;
}

Atom* Residue::getAtom (string atom_name) const{
  //auto it = name_to_atom_map.find(atom_name);
  //if (it == name_to_atom_map.end())
  //  return nullptr;
//  return (*it)->second;
  for(auto const& atom: m_atoms)
    if(atom->getName()==atom_name) return atom;
  return nullptr;
}

std::list<Atom*>& Residue::getAtoms(){
  return m_atoms;
}

const Chain* Residue::getChain() const{
	return m_parentChain;
}

Residue* Residue::getPrevResidue() const {
	return m_prevResidue;
}

Residue* Residue::getNextResidue () const {
	return m_nextResidue;
}

string Residue::getName () const {
	return Name;
}

string Residue::getProperName () const { // should modify this when deal with amber03 force field
	string name = Name;
	
	// Deal with N-/C- terminal residues and not-matching residue names in amber03
	if ( name.size()==4 && (name.substr(0,1)=="N" || ( name.substr(0,1)=="C" && name.substr(1,1)!="Y" )) )
		name = Name.substr(1,3);
	if (name=="HID" || name=="HIP" || name=="HIE") name = "HIS";
	if (name=="HSD" || name=="HSE" || name=="HSP") name = "HIS";
	else if (name=="LYP") name = "LYS";
	else if (name=="CYN" || name=="CYM") name = "CYS";
	else if (name.substr(0,2)=="TE") name = "TES"; // for testing purpose

	// Check names	
	if (name.compare("ALA")==0 || name.compare("VAL")==0 || name.compare("LEU")==0 || name.compare("ILE")==0 
	    || name.compare("PHE")==0 || name.compare("TRP")==0 || name.compare("MET")==0 || name.compare("PRO")==0 
	    || name.compare("ASP")==0 || name.compare("GLU")==0 || name.compare("GLY")==0 || name.compare("SER")==0
	    || name.compare("THR")==0 || name.compare("TYR")==0 || name.compare("CYS")==0 || name.compare("ASN")==0
	    || name.compare("GLN")==0 || name.compare("LYS")==0 || name.compare("ARG")==0 || name.compare("HIS")==0
		|| name.compare("NLE")==0 || name.compare("TES")==0) ;
	else {
		cerr << "Residue " << Id << " has improper name: " << Name << endl;
		exit(1);
	}
	return name;
}

string Residue::getNameType () const {
// Nonpolar: Ala, Val, Leu, Ile, Phe, Trp, Met, Pro
// Polar_acidic: Asp, Glu
// Polar_uncharged: Gly, Ser, Thr, Cys, Tyr, Asn, Gln
// Polar_basic: Lys, Arg, His
	string name = getProperName();
	string type = "Nonpolar";
	if (name.compare("ASP")==0 || name.compare("GLU")==0) 
		type = "Polar_acidic";
	else if (name.compare("GLY")==0 || name.compare("SER")==0 || name.compare("THR")==0 || name.compare("CYS")==0 || name.compare("TYR")==0 || name.compare("ASN")==0 || name.compare("GLN")==0) 
		type = "Polar_uncharged";
	else if (name.compare("LYS")==0 || name.compare("ARG")==0 || name.compare("HIS")==0)
		type = "Polar_basic";
	return type;
}

int Residue::sidechainSize () const {
	int size = 0;
//	for (map<string,Atom*>::const_iterator it= name_to_atom_map.begin(); it != name_to_atom_map.end(); ++it) {
//		if ( (it->second)->isSidechainAtom() )
//			++size;
//	}
  for (auto const& atom: m_atoms)
    if ( atom->isSidechainAtom() ) ++size;
	return size;
}

void Residue::setSSE (string& type) {
	SSE_type = type;
}

string Residue::getSSE () const {
	return SSE_type;
}


bool Residue::isWithinResidueRange( unsigned int resid1, unsigned int resid2 ) const {
	if ( resid2>=resid1 )
		if( getId() >= resid1 && getId() <= resid2 )
			return true;
	return false;
}

bool Residue::isWithinTwoResidueRanges( unsigned int resid1, unsigned int resid2,
					   unsigned int resid3, unsigned int resid4 ) const {
	if ( resid2>=resid1 && resid4>=resid3 )
		if( ( getId() >= resid1 && getId() <= resid2 )
		 || ( getId() >= resid3 && getId() <= resid4 ) )
			return true;
	return false;
}
