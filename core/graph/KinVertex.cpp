#include "KinVertex.h"

#include "Logger.h"

using namespace std;

KinVertex::KinVertex ():
    id(0),
    m_rigidbody(NULL),
    m_parent(NULL),
    isRibose(false)
{
  m_transformation.setIdentity();
}

KinVertex::KinVertex (int id_, Rigidbody* rb_ptr):
    id(id_),
    m_rigidbody(rb_ptr)
{
  m_parent = NULL;
  Visited = false;
  isRibose = false;
  rb_ptr->setVertex(this);
  m_transformation.setIdentity();
}

KinVertex::~KinVertex () {
  for (auto eit=m_edges.begin(); eit!=m_edges.end(); ++eit) {
    delete *eit;
  }
  m_rigidbody->setVertex(NULL);
}

void KinVertex::setParent(KinVertex* v) {
  m_parent = v;
}

void KinVertex::addEdge (unsigned int neighbor_vertex_id, KinEdge *edge) {
  m_edges.push_back( edge );
}

void KinVertex::print() const {
  log() << "KinVertex_" << id << ", id ";
  for (vector<Atom*>::iterator it=m_rigidbody->Atoms.begin(); it!=m_rigidbody->Atoms.end(); ++it)
    log() << (*it)->getId() << "+";
  log() << endl;
}

void KinVertex::forwardPropagate()
{
  transformAtoms();

  for(auto const& edge: m_edges){
    edge->forwardPropagate();
  }
}

void KinVertex::transformAtoms()
{
  //cout<<"KinVertex::transformAtoms - "<<this->id<<endl<<m_transformation<<endl;
  for (auto const& atom: m_rigidbody->Atoms){
    Math3D::Vector3 newPos = m_transformation * atom->m_referencePosition;
    if(atom->getResidue()->getId()==1 && atom->getName()=="CA"){
      cout<<"1/CA trans: "<<endl<<m_transformation<<endl;
      cout<<"1/CA pos:   "<<atom->m_Position<<endl;
      cout<<"1/CA resn:  "<<atom->getResidue()->getName()<<endl;
    }

    //newPos.x += m_transformation.t.x;
    //newPos.y += m_transformation.t.y;
    //newPos.z += m_transformation.t.z;

    atom->m_Position.x = newPos.x;
    atom->m_Position.y = newPos.y;
    atom->m_Position.z = newPos.z;
    if(atom->getResidue()->getId()==1 && atom->getName()=="CA"){
      cout<<"1/CA newpos:"<<atom->m_Position<<endl;
    }
  }
}

void KinVertex::TransformAtomPosition(Math3D::RigidTransform *trsfm){

  for (auto const& atom: m_rigidbody->Atoms){
    Math3D::Vector3 newPos = trsfm->R * atom->m_Position;

    newPos.x += trsfm->t.x;
    newPos.y += trsfm->t.y;
    newPos.z += trsfm->t.z;

    atom->m_Position.x = newPos.x;
    atom->m_Position.y = newPos.y;
    atom->m_Position.z = newPos.z;
  }
}
