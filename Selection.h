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
#ifndef SELECTION_H
#define SELECTION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Residue.h"

using namespace std;

class Selection {
public:
	Selection( );
	Selection( string selection );
	Selection( string selection, string delim );
	void print() const;
	void print( string selName ) const;

	// Mutator and Accessor
	void delim(string delim);
	string delim() const;
	// Mutator and Accessor
	void selection(string selection);
	string selection() const;
	// Mutator and Accessor
	void selectionWords( vector<string> selectionWords );
	vector<string> selectionWords() const;

	vector<Residue*> getSelectedResidues( const Molecule *protein ) const;
	vector<Atom*> getSelectedAtoms( const Molecule *protein );
	vector<Atom*> getSelectedAtoms( const vector<Residue*> residues );

	static vector<string> &split( const string &s, string delim, vector<string> &words );
	static vector<string>  split( const string &s, string delim );
	static vector<int> &split( const string &s, string delim, vector<int> &numbers );
	string &combine( const vector<string> &words, string delim, string &s );
	string  combine( const vector<string> &words, string delim );

private:
	string selection_, delim_;
	vector<string> selectionWords_;
};

#endif // SELECTION_H

