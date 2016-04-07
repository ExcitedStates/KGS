//
// Created by Rasmus Fonseca on 4/7/16.
//

#ifndef KGS_KINTREE_H
#define KGS_KINTREE_H


#include "KinGraph.h"

class RigidbodyTree: public KinGraph {
  // In RigidbodyTree, only edges to m_children are stored in Edges.
  // The bonds associated with edges not necessarily point from smaller Id atom to larger Id atom.
  // Rather, they point from m_parent rigidbody to child rigidbody.
public:
  RigidbodyTree();
  KinVertex *root;
  std::vector< std::pair<Edge*,KinVertex*> > CycleAnchorEdges; // pair<edge,common anchor>; each edge closes a cycle
  int Cycle_DOF_num; // total number of DOFs in the cycles
  int num_DOFs;

  ~RigidbodyTree();
  void print();
  void printForSpringy();
  KinVertex* findCommonAncestor (KinVertex *v1, KinVertex *v2);
};



#endif //KGS_KINTREE_H
