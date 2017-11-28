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

#include "metrics/RMSDnosuper.h"

using namespace std;

namespace metrics {

RMSDnosuper::RMSDnosuper( Selection &selection ) :
    Metric(selection) { }

double RMSDnosuper::distance( Configuration *c1, Configuration *c2 ) {
  const std::vector<Atom*>& atomsRMSD1 = m_selection.getSelectedAtoms(c1->getMolecule());
  const std::vector<Atom*>& atomsRMSD2 = m_selection.getSelectedAtoms(c2->getMolecule());

  if( atomsRMSD1.empty() || atomsRMSD2.empty() ){
    cerr<<"RMSD::distance - Atom-selection given to RMSD metric contained no atoms: "<<m_selection<<endl;
    exit(-1);
  }
  if(atomsRMSD1.size() != atomsRMSD2.size()){
    cerr<<"RMSDnosuper::distance(..)";
    cerr<<" - Configurations have different number of atoms ("<<atomsRMSD1.size()<<" vs "<<atomsRMSD2.size()<<")"<<endl;
    exit(-1);
  }

  vector<Coordinate> p1_atoms;
  c1->updateMolecule();
  for( auto const &a : atomsRMSD1 ) {
    p1_atoms.push_back(a->m_position);
  }

  vector<Coordinate> p2_atoms;
  c2->updateMolecule();
  for( auto const &a : atomsRMSD2 ) {
    p2_atoms.push_back(a->m_position);
  }

  size_t numAtoms = atomsRMSD1.size();
  double sum = 0;
  for( size_t i=0; i < numAtoms; i++ ) {
    sum+=p1_atoms[i].distanceSquared(p2_atoms[i]);
  }

  return sqrt(sum/numAtoms);
}

}
