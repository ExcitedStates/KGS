#ifndef KGS_KINVERTEX_H
#define KGS_KINVERTEX_H


#include "core/graph/RBEdge.h"

class Edge;

class KinVertex {
public:
  const unsigned int id;
  Rigidbody * const Rb_ptr;
  std::vector<Edge*> edges;                ///< Child-edges after spanning tree has been created
  KinVertex *Parent;                       ///< Parent-vertex after spanning tree has been created
  bool Visited;   ///< When finding common ancestor, vertices are marked as visited up to the root
  Math3D::RigidTransform transformation;   ///< The rigid body transformation matrix of this rigid group

  bool isRibose;

  KinVertex();
  KinVertex(int id, Rigidbody* rb);
  virtual ~KinVertex();

  void addEdge(unsigned int neighbor_vertex_id, Edge *edge);
  Edge* findEdge(KinVertex* v) const;
  virtual void setParent(KinVertex* v);
  void print();
  void TransformAtomPosition(Math3D::RigidTransform *trsfm);
};


#endif //KGS_KINVERTEX_H
