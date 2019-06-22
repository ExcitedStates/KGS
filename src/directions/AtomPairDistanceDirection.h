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

#ifndef KGS_ATOMPAIRDISTANCEDIRECTION_H
#define KGS_ATOMPAIRDISTANCEDIRECTION_H


#include <gsl/gsl_vector.h>

#include <tuple>
#include <vector>

#include "core/graph/KinTree.h"
#include "directions/Direction.h"
#include "Selection.h"

/**
 * A direction where the gradient seeks to drive configurations toward a desired distance between a particular pair of
 * atoms. This class is conceptually similar to `RelativeMSDDirection` and to a lesser degree `LSNrelativeDirection`,
 * but only takes a single pair of atoms and desired distance.
 */
class AtomPairDistanceDirection: public Direction {
 public:
  AtomPairDistanceDirection(const std::tuple<Atom*, Atom*, double> &relativeDistances);

 protected:
  void computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret);

 private:
  const std::tuple<Atom*, Atom*, double> &m_relativeDistance;
};


#endif //KGS_ATOMPAIRDISTANCEDIRECTION_H
