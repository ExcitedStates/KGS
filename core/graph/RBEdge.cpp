
#include "RBEdge.h"
#include "Logger.h"

using namespace std;

Edge::Edge(KinVertex *startv, KinVertex *endv, Bond *bond) :
    StartVertex(startv),
    EndVertex(endv),
    m_bond(bond)
{
  DOF_id = -1;
  Cycle_DOF_id = -1;
}

void Edge::print () const {
  log() << "(" << StartVertex->id << "," << EndVertex->id;
  log() << ",[" << m_bond->Atom1->getResidue()->getId()<<m_bond->Atom1->getName() << ",";
  log() << m_bond->Atom2->getResidue()->getId() << m_bond->Atom2->getName() << "]";
  log() << ",[" << m_bond->Atom1->getId() << "," << m_bond->Atom2->getId() << ",";
  log() << (m_bond->BondType=="HB"?"HB":"CV") << "]";
  log() << ",DOF_id=" << DOF_id << ",Cycle_DOF_id=" << Cycle_DOF_id << ")" << endl;
}

void Edge::printVerbose() const{
  log() << "\t Edge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t StartVertex: " << (StartVertex->Rb_ptr->isMainchainRb()?"isOnMainchain":"isNOTonMainchain") << endl;
  //log() << "\t\t Direction: Res1_ID Atom1_Name Atom1_ID ---------> Res2_ID Atom2_Name Atom2_ID " << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "\t\t Bond type: " << (m_bond->BondType=="HB"?"HB":"CV") << ", " << endl;
  log() << "\t\t DOF_id= " << DOF_id << "  Cycle_DOF_id= " << Cycle_DOF_id << endl;
  log() << "---------------------------" << endl;
}

void Edge::printShort() const{
  log() << "\t Edge: " << endl;
  log() << "\t\t StartVertex: " << StartVertex->id << "  EndVertex: " << EndVertex->id << endl;
  log() << "\t\t Direction: " << m_bond->Atom1->getResidue()->getId() << " " << m_bond->Atom1->getId() << " " << m_bond->Atom1->getName() <<
  "\t ---------> \t"
  << m_bond->Atom2->getResidue()->getId() << " " << m_bond->Atom2->getId() << " " << m_bond->Atom2->getName() << endl;
  log() << "---------------------------" << endl;
}

void Edge::printHTML() const {
  int label = EndVertex->Rb_ptr->size();
  //int label = EndVertex->id;
  log() << "var vertex_" << EndVertex->id << " = graph.newNode({label: \"" << label << "\", fill: \"#000000\"});" << endl;
  log() << "graph.newEdge(vertex_" << StartVertex->id << ", vertex_" << EndVertex->id << ");" << endl;
}

void Edge::printHTMLRoot() const{
  log() << "var vertex_" << StartVertex->id << " = graph.newNode({label: \"*"
  << StartVertex->id << "*\", fill: \"#EE2222\"});" << endl;
}

Bond *Edge::getBond() const {
  return m_bond;
}



ostream& operator<<(ostream& os, const Edge& e){
  os<<"Edge["<<e.getBond()->Atom1->getName()<<", "<<e.getBond()->Atom2->getName()<<"]";
  return os;
}

ostream& operator<<(ostream& os, const Edge* e){
  os<<*e;
  return os;
}
