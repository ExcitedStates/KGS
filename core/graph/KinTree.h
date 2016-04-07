#ifndef KGS_KINTREE_H
#define KGS_KINTREE_H


#include <set>
#include "KinGraph.h"

/**
 * In KinTree, only edges to children are stored in `m_edges`.
 * The bonds associated with m_edges not necessarily point from smaller Id atom to larger Id atom.
 * Rather, they point from m_parent rigidbody to child rigidbody.
 */
class KinTree: public KinGraph
{
 public:
  KinVertex *root;
  std::vector< std::pair<KinEdge*,KinVertex*> > CycleAnchorEdges; // pair<edge,common anchor>; each edge closes a cycle
  int m_numCycleDOFs; // total number of DOFs in the cycles
  int m_numDOFs;
  std::vector<DOF*> m_dofs;

  KinTree();
  ~KinTree();
  void print();
  void printForSpringy();
  KinVertex* findCommonAncestor (KinVertex *v1, KinVertex *v2);
  void collectDOFs();
private:
  void collectDOFs(KinVertex* v);


};



#endif //KGS_KINTREE_H
