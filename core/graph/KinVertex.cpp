#include "KinVertex.h"

#include "Logger.h"

using namespace std;

KinVertex::KinVertex ():
    id(0),
    Rb_ptr(NULL)
{
  Parent = NULL;
  isRibose = false;
}

KinVertex::KinVertex (int id_, Rigidbody* rb_ptr):
    id(id_),
    Rb_ptr(rb_ptr)
{
  Parent = NULL;
  Visited = false;
  isRibose = false;
  rb_ptr->setVertex(this);
}

KinVertex::~KinVertex () {
  for (auto eit=edges.begin(); eit!=edges.end(); ++eit) {
    delete *eit;
  }
  Rb_ptr->setVertex(NULL);
}

void KinVertex::setParent(KinVertex* v) {
  Parent = v;
}

void KinVertex::addEdge (unsigned int neighbor_vertex_id, Edge *edge) {
  edges.push_back( edge );
}

void KinVertex::print () {
  log() << "KinVertex_" << id << ", id ";
  for (vector<Atom*>::iterator it=Rb_ptr->Atoms.begin(); it!=Rb_ptr->Atoms.end(); ++it)
    log() << (*it)->getId() << "+";
  log() << endl;
}

void KinVertex::TransformAtomPosition(Math3D::RigidTransform *trsfm){

  for (auto const& atom: Rb_ptr->Atoms){
    Math3D::Vector3 newPos = trsfm->R * atom->m_Position;

    newPos.x += trsfm->t.x;
    newPos.y += trsfm->t.y;
    newPos.z += trsfm->t.z;

    atom->m_Position.x = newPos.x;
    atom->m_Position.y = newPos.y;
    atom->m_Position.z = newPos.z;
  }
}
