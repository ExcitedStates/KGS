#include <core/dofs/TorsionDOF.h>
#include <core/dofs/GlobalRotateDOF.h>
#include <core/dofs/GlobalTranslateDOF.h>
#include "KinEdge.h"
#include "Logger.h"

using namespace std;

KinEdge::KinEdge(KinVertex *startv, KinVertex *endv, Bond *bond, int dof_id) :
    StartVertex(startv),
    EndVertex(endv),
    m_bond(bond),
    m_dof(NULL),
    DOF_id(dof_id)
{
  Cycle_DOF_id = -1;
}

DOF* KinEdge::setDOF(DOF* dof) {
  if(m_dof){
    std::cerr<<
  }

}

void KinEdge::print () const {
  log() << "(" << StartVertex->id << "," << EndVertex->id;
  log() << ",[" << m_bond->Atom1->getResidue()->getId()<<m_bond->Atom1->getName() << ",";
  log() << m_bond->Atom2->getResidue()->getId() << m_bond->Atom2->getName() << "]";
  log() << ",[" << m_bond->Atom1->getId() << "," << m_bond->Atom2->getId() << ",";
  log() << (m_bond->BondType=="HB"?"HB":"CV") << "]";
  log() << ",DOF_id=" << DOF_id << ",Cycle_DOF_id=" << Cycle_DOF_id << ")" << endl;
}

void KinEdge::printVerbose() const{
  log() << "\t KinEdge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t StartVertex: " << (StartVertex->m_rigidbody->isMainchainRb()?"isOnMainchain":"isNOTonMainchain") << endl;
  //log() << "\t\t Direction: Res1_ID Atom1_Name Atom1_ID ---------> Res2_ID Atom2_Name Atom2_ID " << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "\t\t Bond type: " << (m_bond->BondType=="HB"?"HB":"CV") << ", " << endl;
  log() << "\t\t DOF_id= " << DOF_id << "  Cycle_DOF_id= " << Cycle_DOF_id << endl;
  log() << "---------------------------" << endl;
}

void KinEdge::printShort() const{
  log() << "\t KinEdge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "---------------------------" << endl;
}

void KinEdge::printHTML() const {
  int label = EndVertex->m_rigidbody->size();
  //int label = EndVertex->id;
  log() << "var vertex_" << EndVertex->id << " = graph.newNode({label: \"" << label << "\", fill: \"#000000\"});" << endl;
  log() << "graph.newEdge(vertex_" << StartVertex->id << ", vertex_" << EndVertex->id << ");" << endl;
}

void KinEdge::printHTMLRoot() const{
  log() << "var vertex_" << StartVertex->id << " = graph.newNode({label: \"*"
  << StartVertex->id << "*\", fill: \"#EE2222\"});" << endl;
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
