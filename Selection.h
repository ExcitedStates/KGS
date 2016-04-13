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


class Selection {
public:
	Selection( );
	Selection( std::string selection );
	Selection( std::string selection, std::string delim );
	void print() const;
	void print( std::string selName ) const;

	// Mutator and Accessor
	void delim(std::string delim);
	std::string delim() const;
	// Mutator and Accessor
	void selection(std::string selection);
	std::string selection() const;
	// Mutator and Accessor
	void selectionWords( std::vector<std::string> selectionWords );
	std::vector<std::string> selectionWords() const;

	std::vector<Residue*> getSelectedResidues( const Molecule *protein ) const;
	std::vector<Atom*> getSelectedAtoms( const Molecule *protein );
	std::vector<Atom*> getSelectedAtoms( const std::vector<Residue*> residues );

	static std::vector<std::string> &split( const std::string &s, std::string delim, std::vector<std::string> &words );
	static std::vector<std::string>  split( const std::string &s, std::string delim );
	static std::vector<int> &split( const std::string &s, std::string delim, std::vector<int> &numbers );
	std::string &combine( const std::vector<std::string> &words, std::string delim, std::string &s );
	std::string  combine( const std::vector<std::string> &words, std::string delim );

private:
	std::string selection_, delim_;
	std::vector<std::string> selectionWords_;
};

#endif // SELECTION_H

