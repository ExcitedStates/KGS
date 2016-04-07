/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

        Permission is hereby granted, free of charge, to any person obtaining a copy of
        this software and associated documentation files (the "Software"), to deal in
        the Software without restriction, including without limitation the rights to
        use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
        of the Software, and to permit persons to whom the Software is furnished to do
        so, subject to the following conditions:

        This entire text, including the above copyright notice and this permission notice
        shall be included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
        OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.


*/
#ifndef RIGIDBODYGRAPH_H
#define RIGIDBODYGRAPH_H

#include <vector>
#include <map>

#include "math/MathUtility.h"
#include "core/Bond.h"
#include "Rigidbody.h"
#include "math3d/primitives.h"

class RigidbodyGraphVertex;
class RigidbodyTree;

class Edge {
 public:
  RigidbodyGraphVertex *StartVertex;
  RigidbodyGraphVertex *EndVertex;
  //bool Hbond;
  bool Shared;
  unsigned int RigidbodyGraphCycleCluster_id;
  int DOF_id; // Start from 0. If the edge is not a DOF, its DOF_id is -1.
  int Cycle_DOF_id; // IDs of DOFs in cycles only. Start from 0. If the edge is not a cycle dof, the value is -1.

  Edge(RigidbodyGraphVertex *startv, RigidbodyGraphVertex *endv, Bond * m_bond);
  void print();
  void printVerbose();
  void printShort();
  void printHTML();
  void printHTMLRoot();

  Bond *getBond() const;

 private:
  Bond * const m_bond;
};

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

class RigidbodyGraphCycle {
 public:
  std::vector<Edge*> CycleEdges;
  unsigned int IndependentLength; // number of independent DOFs in the cycle; H-bond is NOT a DOF!

  RigidbodyGraphCycle();
  ~RigidbodyGraphCycle() {};
};

class RigidbodyGraphCycleCluster {
 public:
  std::vector< RigidbodyGraphCycle > Cycles;

  void print();
};

class RigidbodyGraph {
 public:
  std::map<unsigned int, RigidbodyGraphVertex*> Vertex_map;
  std::list< std::pair< unsigned int, RigidbodyGraphVertex*> > m_sortedVertices;
  std::vector<Edge*> Edges;
  std::map<unsigned int, RigidbodyGraphCycleCluster> CycleClusters; // cluster_id -> cluster
  unsigned int MaxCycleClusterId;

  RigidbodyGraph ();
  ~RigidbodyGraph ();
  RigidbodyGraphVertex* addVertex (unsigned int vertex_id, Rigidbody *rb, bool flexibleSugar);
  RigidbodyGraphVertex* getVertex (int rb_id);
  void addEdge (RigidbodyGraphVertex *vertex1, RigidbodyGraphVertex *vertex2, Bond * bond);
  Edge* addEdgeDirected (RigidbodyGraphVertex *vertex1, RigidbodyGraphVertex *vertex2, Bond * bond, int DOF_id); // Add a directed edge from rb_id1 to rb_id2

  bool hasVertex (int rb_id);
  void print ();
  void findCycleClusters();
};

class RigidbodyTree: public RigidbodyGraph {
  // In RigidbodyTree, only edges to m_children are stored in Edges.
  // The bonds associated with edges not necessarily point from smaller Id atom to larger Id atom.
  // Rather, they point from m_parent rigidbody to child rigidbody.
 public:
  RigidbodyTree();
  RigidbodyGraphVertex *root;
  std::vector< std::pair<Edge*,RigidbodyGraphVertex*> > CycleAnchorEdges; // pair<edge,common anchor>; each edge closes a cycle
  int Cycle_DOF_num; // total number of DOFs in the cycles
  int num_DOFs;

  ~RigidbodyTree();
  void print();
  void printForSpringy();
  RigidbodyGraphVertex* findCommonAncestor (RigidbodyGraphVertex *v1, RigidbodyGraphVertex *v2);
};

std::ostream& operator<<(std::ostream& os, const Edge& e);
std::ostream& operator<<(std::ostream& os, const Edge* e);

#endif
