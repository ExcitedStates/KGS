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

void KinEdge::printVerbose() const{
  log() << "\t KinEdge: " << endl;
  log() << "\t\t StartVertex: " << (StartVertex->m_rigidbody->isMainchainRb()?"isOnMainchain":"isNOTonMainchain") << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "\t\t Bond type: " << (m_bond->BondType=="HB"?"HB":"CV") << ", " << endl;
  log() << "---------------------------" << endl;
}

void KinEdge::printShort() const{
  log() << "\t KinEdge: " << endl;
//  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "---------------------------" << endl;
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


ostream& operator<<(ostream& os, const KinEdge& e){
  os<<"KinEdge["<<e.getBond()->Atom1->getName()<<", "<<e.getBond()->Atom2->getName()<<"]";
  return os;
}

ostream& operator<<(ostream& os, const KinEdge* e){
  os<<*e;
  return os;
}
