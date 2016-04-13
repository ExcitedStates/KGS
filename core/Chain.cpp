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

#include "Chain.h"
#include "Residue.h"
#include "Molecule.h"
#include "Atom.h"

using namespace std;

Chain::Chain (): Name("UNKNOWN")
{
}

Chain::Chain(const string& name, Molecule * parent_protein):
		Name(name)
{
	Parent_protein = parent_protein;
}

Chain::~Chain() {
	// delete all residues
	for (vector<Residue*>::iterator it=Residue_list.begin(); it!=Residue_list.end(); ++it) {
		delete (*it);
	}
}

const string& Chain::getName () const {
	return Name;
}

Atom* Chain::addAtom (const std::string& resName,
                     const int& resId,
                     const std::string& atomName,
                     const int& atomId,
                     const Coordinate& position )
{
	Residue* res = getResidue(resId);
	if (res == nullptr) { // this is a new residue
		res = addResidue(resName,resId);
	}
	return res->addAtom(atomName, atomId, position);
}

std::vector<Residue*>& Chain::getResidues()
{
  return Residue_list;
}

Residue* Chain::getResidue (int res_id){
	for (vector<Residue*>::iterator it=Residue_list.begin(); it!=Residue_list.end(); ++it) {
		if ((*it)->getId() == res_id) {
			return (*it);
		}
	}
	return nullptr;
}

Residue* Chain::addResidue (const string& resName, const int& resId) {
	Residue* res = new Residue(resName,resId,this);
	Residue_list.push_back(res);
	if (Residue_list.size()>1) { // there are at least 2 residues
		Residue* previous_res = Residue_list.at(Residue_list.size()-2);
		previous_res->setNextResidue(res);
	}
	return res;
}

Molecule * Chain::getProtein() const {
	return Parent_protein;
}

void Chain::printSummaryInfo() const {
	cout << "Chain " << Name << endl;
	for (vector<Residue*>::const_iterator it=Residue_list.begin(); it!=Residue_list.end(); ++it) {
		(*it)->printSummaryInfo();
	}
}

void Chain::print() const {
	cout << "Chain " << Name << endl;
	for (vector<Residue*>::const_iterator it=Residue_list.begin(); it!=Residue_list.end(); ++it) {
		cout<<"\t";
		(*it)->print();
	}
}

int Chain::size () const {
	return Residue_list.size();
}

