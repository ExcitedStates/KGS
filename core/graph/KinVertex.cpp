#include "KinVertex.h"

#include "Logger.h"

using namespace std;

KinVertex::KinVertex (Rigidbody* rb_ptr):
    m_rigidbody(rb_ptr)
{
  m_parent = nullptr;
  Visited = false;

  if(rb_ptr!=nullptr)
    rb_ptr->setVertex(this);

  m_transformation.setIdentity();
}

KinVertex::~KinVertex () {
  for (auto eit=m_edges.begin(); eit!=m_edges.end(); ++eit) {
    delete *eit;
  }
  if(m_rigidbody)
    m_rigidbody->setVertex(nullptr);
}

void KinVertex::setParent(KinVertex* v) {
  m_parent = v;
}

void KinVertex::addEdge (KinEdge *edge) {
  m_edges.push_back( edge );
}

KinEdge* KinVertex::findEdge(const KinVertex* v) const
{
  for(auto const& edge: m_edges){
    if( edge->EndVertex==v )
      return edge;
  }
  return nullptr;
}


void KinVertex::print() const {
  log() << "KinVertex";
  if(m_rigidbody)
    for (vector<Atom*>::iterator it=m_rigidbody->Atoms.begin(); it!=m_rigidbody->Atoms.end(); ++it)
      log() << (*it)->getId() << "+";
  log() << endl;
}

void KinVertex::forwardPropagate()
{
  for(auto const& edge: m_edges){
    edge->forwardPropagate();
  }

  //Apply transformation AFTER propagation. This is important.
  transformAtoms();
}

void KinVertex::transformAtoms()
{
  if(m_rigidbody==nullptr) return;

  for (auto const& atom: m_rigidbody->Atoms){
    Math3D::Vector3 newPos = m_transformation * atom->m_referencePosition;

    atom->m_Position.x = newPos.x;
    atom->m_Position.y = newPos.y;
    atom->m_Position.z = newPos.z;
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
