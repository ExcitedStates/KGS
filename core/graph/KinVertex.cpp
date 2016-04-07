#include <core/Rigidbody.h>
#include "KinVertex.h"
#include "RBEdge.h"

class RigidbodyGraphVertex {
 public:
  const unsigned int id;
  Rigidbody * const Rb_ptr;
  std::vector<Edge*> edges;     ///< Child-edges after spanning tree has been created
  RigidbodyGraphVertex *Parent; ///< Parent-vertex after spanning tree has been created
  bool Visited;                 ///< When finding common ancestor, vertices are marked as visited up to the root
  Math3D::RigidTransform transformation;  ///< The rigid body transformation matrix of this rigid group

  bool isRibose;

  RigidbodyGraphVertex();
  RigidbodyGraphVertex(int id, Rigidbody* rb);
  virtual ~RigidbodyGraphVertex();

  void addEdge(unsigned int neighbor_vertex_id, Edge *edge);
  Edge* findEdge(RigidbodyGraphVertex* v) const;
  virtual void setParent(RigidbodyGraphVertex* v);
  void print();
  void TransformAtomPosition(Math3D::RigidTransform *trsfm);
};
