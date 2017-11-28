/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

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
