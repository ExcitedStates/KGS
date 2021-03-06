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


#ifndef KGS_MSDDIRECTION_H
#define KGS_MSDDIRECTION_H


#include <core/graph/KinTree.h>
#include "Direction.h"
#include "Selection.h"
#include "metrics/RMSD.h"

class MSDDirection: public Direction {
 public:
  MSDDirection(Selection& resNetwork, bool alignAlways = false);

 protected:
  void computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret);

 private:
  /** Build the `m_sortedVertices`. */
  void collectVerticesPostorder(KinVertex*);

  /** A preprocessed list of all vertices (except the m_root) sorted according to a post-order traversal of the tree. */
  std::list< KinVertex* > m_sortedVertices;

  KinTree* m_preprocessedTree;

  Selection& m_resNetwork;
  metrics::RMSD m_rmsd;
  bool m_alignAlways;
};


#endif //KGS_MSDDIRECTION_H
