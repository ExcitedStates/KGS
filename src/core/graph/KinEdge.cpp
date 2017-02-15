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
    log() << "KinEdge[" << m_bond->Atom1->getResidue()->getId() << m_bond->Atom1->getName() << ",";
    log() << m_bond->Atom2->getResidue()->getId() << m_bond->Atom2->getName()<<",";
    log() << "ID" << m_bond->Atom1->getId() << ",ID" << m_bond->Atom2->getId() << ",";
    log() << (m_bond->BondType == "HB" ? "HB" : "CV") << "]" << endl;
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
  minID1 = edge1->m_bond->Atom1->getId();
  if (edge1->m_bond->Atom2->getId() < minID1 ){
    maxID1 = minID1;
    minID1 = edge1->m_bond->Atom2->getId();
  }
  else
    maxID1 = edge1->m_bond->Atom2->getId();

  minID2 = edge2->m_bond->Atom1->getId();
  if (edge2->m_bond->Atom2->getId() < minID2 ){
    maxID2 = minID2;
    minID2 = edge2->m_bond->Atom2->getId();
  }
  else
    maxID2 = edge2->m_bond->Atom2->getId();
  //Sort
  if(minID1 < minID2 )
    return true;
  if(minID1 > minID2)
    return false;
  return maxID1 < maxID2;
}


ostream& operator<<(ostream& os, const KinEdge& e){
  //os<<"KinEdge["<<e.getBond()->Atom1->getName()<<", "<<e.getBond()->Atom2->getName()<<"]";
  os<<"KinEdge["<<e.getBond()->Atom1<<":"<<e.getBond()->Atom1->getId()<<", "<<e.getBond()->Atom2<<":"<<e.getBond()->Atom2->getId()<<"]";
  return os;
}

ostream& operator<<(ostream& os, const KinEdge* e){
  os<<*e;
  return os;
}
