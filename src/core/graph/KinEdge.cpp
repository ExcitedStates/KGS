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

#include <core/dofs/TorsionDOF.h>
#include <core/dofs/GlobalRotateDOF.h>
#include <core/dofs/GlobalTranslateDOF.h>
#include "KinEdge.h"
#include "Logger.h"

using namespace std;

KinEdge::KinEdge(KinVertex *startv, KinVertex *endv, Bond *bond) :
    StartVertex(startv),
    EndVertex(endv),
    m_bond(bond),
    m_dof(nullptr)
//    DOF_id(dof_id)
{
//  Cycle_DOF_id = -1;
}
KinEdge::~KinEdge(){
  if(m_dof){
    delete m_dof;
  }
}

void KinEdge::setDOF(DOF* dof) {
  if(m_dof){
    std::cerr<<"KinEdge::setDOF - m_dof already set."<<std::endl;
    exit(-1);
  }

   m_dof = dof;
}

void KinEdge::print () const {
  if(m_bond==nullptr)
    log()<<"KinEdge[global]"<<endl;
  else {
    log() << "KinEdge[" << m_bond->m_atom1->getResidue()->getId() << m_bond->m_atom1->getName() << ",";
    log() << m_bond->m_atom2->getResidue()->getId() << m_bond->m_atom2->getName()<<",";
    log() << "ID" << m_bond->m_atom1->getId() << ",ID" << m_bond->m_atom2->getId() << ",";
    log() << (m_bond->m_bondType == "HB" ? "HB" : "CV") << "]" << endl;
  }
}

Bond *KinEdge::getBond() const {
  return m_bond;
}

DOF *KinEdge::getDOF() const {
  return m_dof;
}

void KinEdge::forwardPropagate()
{
  m_dof->updateEndVertexTransformation();
  EndVertex->forwardPropagate();
}

///Compare IDs of two bonds, used to sort them, lowest ID goes first
bool KinEdge::compareIDs(KinEdge* edge1, KinEdge* edge2) {
  //Determine min/max IDs in case bonds were not from min to max ID (happens for hbonds)
  int minID1, maxID1, minID2, maxID2;
  minID1 = edge1->m_bond->m_atom1->getId();
  if (edge1->m_bond->m_atom2->getId() < minID1 ){
    maxID1 = minID1;
    minID1 = edge1->m_bond->m_atom2->getId();
  }
  else
    maxID1 = edge1->m_bond->m_atom2->getId();

  minID2 = edge2->m_bond->m_atom1->getId();
  if (edge2->m_bond->m_atom2->getId() < minID2 ){
    maxID2 = minID2;
    minID2 = edge2->m_bond->m_atom2->getId();
  }
  else
    maxID2 = edge2->m_bond->m_atom2->getId();
  //Sort
  if(minID1 < minID2 )
    return true;
  if(minID1 > minID2)
    return false;
  return maxID1 < maxID2;
}


ostream& operator<<(ostream& os, const KinEdge& e){
  //os<<"KinEdge["<<e.getBond()->Atom1->getName()<<", "<<e.getBond()->m_atom2->getName()<<"]";
  os<<"KinEdge["<<e.getBond()->m_atom1<<":"<<e.getBond()->m_atom1->getId()<<", "<<e.getBond()->m_atom2<<":"<<e.getBond()->m_atom2->getId()<<"]";
  return os;
}

ostream& operator<<(ostream& os, const KinEdge* e){
  os<<*e;
  return os;
}
