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
  if(rb_ptr!=NULL)
    rb_ptr->setVertex(this);
  m_transformation.setIdentity();
}

KinVertex::~KinVertex () {
  for (auto eit=m_edges.begin(); eit!=m_edges.end(); ++eit) {
    delete *eit;
  }
  if(m_rigidbody)
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
  if(m_rigidbody)
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
  if(!m_rigidbody) return;

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
