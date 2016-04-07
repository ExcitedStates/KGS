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

#ifndef KGS_RBEDGE_H
#define KGS_RBEDGE_H

#include "KinGraph.h"

/**
 * An edge connecting two vertices in the kinematic tree.
 *
 * An edge is responsible for holding its two end-vertices and reference to a DOF object (though one DOF
 * can be associated with more Edges - see documentation of PentamerDOF). The edge is not responsible for
 * freeing any of these.
 */
class Edge {
 public:
  Edge(KinVertex * startv, KinVertex * endv, Bond * m_bond);

  KinVertex *StartVertex;
  KinVertex *EndVertex;

  int DOF_id; // Start from 0. If the edge is not a DOF, its DOF_id is -1.
  int Cycle_DOF_id; // IDs of DOFs in cycles only. Start from 0. If the edge is not a cycle dof, the value is -1.


  void print() const;

  void printVerbose() const;
  void printShort() const;
  void printHTML() const;
  void printHTMLRoot() const;

  Bond *getBond() const;

 private:
  Bond * const m_bond;
};

std::ostream& operator<<(std::ostream& os, const Edge& e);
std::ostream& operator<<(std::ostream& os, const Edge* e);


#endif //KGS_RBEDGE_H
