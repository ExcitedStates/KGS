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
//  int m_numCycleDOFs; // total number of DOFs in the cycles

  KinTree();
  ~KinTree();

  size_t getNumDOFs() const;

  size_t getNumCycleDOFs() const;

  DOF* getDOF(unsigned int idx) const;

  DOF* getCycleDOF(unsigned int idx) const;

  void addCycleDOF(DOF* dof);

  void print() const;

  KinVertex* findCommonAncestor (KinVertex *v1, KinVertex *v2);

  /** Collects all DOFs and indexes them. This is performed after the entire tree is built.  */
  void collectDOFs();

 private:
  void collectDOFs(KinVertex* v);

  std::vector<DOF*> m_cycleDOFs;
  std::vector<DOF*> m_dofs;

  /**
   * A preprocessed list of vertices (and rigidbody id's (?)) sorted according to a post-order traversal
   * of the tree. Used by `MSDDirection::computeGradient`.
   * TODO: Move to MSDDirection.
   */
//  std::list< std::pair< unsigned int, KinVertex*> > m_sortedVertices;

};



#endif //KGS_KINTREE_H
