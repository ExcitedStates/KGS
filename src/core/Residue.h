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

#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <map>
#include <list>
#include "core/Coordinate.h"


class Atom;
class Chain;

class Residue {
 public:
  Residue();
  Residue(const std::string& name, const int& id, const Chain* parent_chain, const std::string& sse_type="UNKNOWN");
  ~Residue();
  bool isEmpty ();
  void setLastResidue(Residue* last);
  void setNextResidue(Residue* next);
  int getId () const;
  Atom* addAtom (
      const bool& hetatm,
      const std::string& atomName,
      const int& atomId,
      const Coordinate& position);
  Atom* getAtom (std::string atom_name) const;
  const std::list<Atom*>& getAtoms() const;
  const Chain* getChain () const;
  void printSummaryInfo() const;
  void print() const;
  Residue* getPrevResidue() const;
  Residue* getNextResidue() const;
  std::string getName () const;
  std::string getProperName () const;
  std::string getNameType () const;
  int sidechainSize () const;
  void setSSE (std::string& type);
  std::string getSSE () const;


  bool isWithinResidueRange( unsigned int resid1, unsigned int resid2 ) const;
  bool isWithinTwoResidueRanges( unsigned int resid1, unsigned int resid2, unsigned int resid3, unsigned int resid4 ) const;

 private:
  const std::string Name;
  const int Id;
  const Chain* m_parentChain;
  std::string SSE_type;
  Residue* m_prevResidue;
  Residue* m_nextResidue;
  //std::map<std::string,Atom*> name_to_atom_map;
  std::list<Atom*> m_atoms;


};

#endif
