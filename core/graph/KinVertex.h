#ifndef KGS_KINVERTEX_H
#define KGS_KINVERTEX_H


#include "core/graph/KinEdge.h"

class KinEdge;

class KinVertex {
public:
  Rigidbody * const m_rigidbody;
  std::vector<KinEdge*> m_edges;             ///< Child-m_edges after spanning tree has been created
  KinVertex *m_parent;                       ///< m_parent-vertex after spanning tree has been created
  //TODO: Visited should be a local variable, not accessible to everyone here
  bool Visited;   ///< When finding common ancestor, vertices are marked as visited up to the m_root
  Math3D::RigidTransform m_transformation;   ///< The transformation to apply to atoms in the rigid body

  KinVertex(Rigidbody* rb=nullptr);
  virtual ~KinVertex();

  void addEdge(KinEdge *edge);
  KinEdge* findEdge(const KinVertex* v) const;
  virtual void setParent(KinVertex* v);
  void print() const;

  void forwardPropagate();
private:
  void transformAtoms();
};


#endif //KGS_KINVERTEX_H
