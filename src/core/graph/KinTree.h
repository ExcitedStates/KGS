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
#ifndef KGS_KINTREE_H
#define KGS_KINTREE_H


#include <set>
#include "KinGraph.h"

/**
 * A kinematic tree (or more correctly set of trees / forest) that defines forward kinematics in the graph. Edges in the
 * graph that are not part of the tree impose constraints.
 * The bonds associated with edges point from parent rigidbody to child rigidbody.
 */
class KinTree: public KinGraph
{
 public:
  KinVertex *m_root;
  std::vector< std::pair<KinEdge*,KinVertex*> > m_cycleAnchorEdges; // pair<edge,common anchor>; each edge closes a cycle

  /** Build a tree that spans the rigid bodies with the specified roots */
  KinTree( const std::vector<Rigidbody*>& rigidbodies, const std::vector<Atom*>& roots );

  ~KinTree();

  /** Return total number of DOFs, both in cycles and free */
  size_t getNumDOFs() const;

  /** Return number of DOFs in cycles only */
  size_t getNumCycleDOFs() const;

  /** Return the DOF with the specified index */
  DOF* getDOF(unsigned int idx) const;

  /** Return the DOF with the specified index */
  DOF* getCycleDOF(unsigned int idx) const;

  void print() const;


  /** Add a directed edge from `vertex1` to `vertex2` */
  KinEdge* addEdgeDirected (KinVertex *vertex1, KinVertex *vertex2, Bond * bond);

  /** Find the lowest common ancestor for v1 and v2 */
  KinVertex* findCommonAncestor (KinVertex *v1, KinVertex *v2);

 private:

  void addCycleDOF(DOF* dof);

  /** Collects all DOFs and indexes them. This is performed after the entire tree is built.  */
  void collectDOFs();

  void collectDOFs(KinVertex* v);

  std::vector<DOF*> m_cycleDOFs;
  std::vector<DOF*> m_dofs;

};


#endif //KGS_KINTREE_H
