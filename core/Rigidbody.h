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
#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include <vector>

#include "core/Atom.h"
#include "core/Bond.h"

class KinVertex;

class Rigidbody {
 public:
  //unsigned int RbId;
  std::vector<Atom*> Atoms;
  std::vector<Bond *> Bonds;

  Rigidbody();
  Rigidbody(unsigned int id);
  ~Rigidbody();

  void id(unsigned int id);
  unsigned int id() const;
  void setMainchainRb();
  bool isMainchainRb() const;

  bool isWithinResidueRange( unsigned int resid1, unsigned int resid2 ) const;
  bool isWithinTwoResidueRanges( unsigned int resid1, unsigned int resid2, unsigned int resid3, unsigned int resid4 ) const;

  Atom* getAtom(std::string name);
  void addAtom(Atom* atom);
  void addBond(Bond * bond);
  void makeBiggerRigidBody( Rigidbody* rb);
  void makeBiggerRigidBody( Rigidbody* rb, Bond * bond);
  int size() const;
  void print() const;
  void printAtoms() const;
  void printAtomsBonds() const;
  bool containsResidue( Residue* res ) const;
  bool containsAtom (Atom* atom) const;
  //bool containsAtomAtPosition( const clipper::Coord_orth& pos ) const;
  bool containsMainchainAtoms() const;
  bool containsAlongMainchainAtoms() const;

  void setVertex (KinVertex* vertex);
  KinVertex* getVertex();

 private:
  unsigned int rbId_;
  bool isMainchainRb_;
  KinVertex* m_rbVertex;
};

std::ostream& operator<<(std::ostream& os, const Rigidbody& rb);
std::ostream& operator<<(std::ostream& os, const Rigidbody* rb);

#endif
