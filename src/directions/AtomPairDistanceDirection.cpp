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
#include "AtomPairDistanceDirection.h"

#include <cmath>
#include <stack>

#include "metrics/RMSD.h"
#include "Logger.h"


using std::vector;
using std::tuple;

AtomPairDistanceDirection::AtomPairDistanceDirection(const std::tuple<Atom*, Atom*, double> &relativeDistance):
    m_relativeDistance(relativeDistance)
{
}

void AtomPairDistanceDirection::computeGradient(Configuration* conf, Configuration* c_target, gsl_vector* ret) {
  assert(c_target == nullptr);
  assert(conf->getMolecule()->m_spanningTree->getNumDOFs() == ret->size);

  Molecule *mol = conf->updatedMolecule();

  vector<double> gradient(ret->size, 0.0);
  int dofCount = 0;

  Atom* a1 = std::get<0>(m_relativeDistance);
  Atom* a2 = std::get<1>(m_relativeDistance);
  double desiredDist = std::get<2>(m_relativeDistance);

  Math3D::Vector3 diff = a2->m_position - a1->m_position;
  double dist = diff.norm();
  diff *= (dist - desiredDist) / dist;

  KinVertex* v1 = a1->getRigidbody()->getVertex();
  KinVertex* v2 = a2->getRigidbody()->getVertex();
  assert(v1 != v2);
  KinVertex* nca = mol->m_spanningTree->findCommonAncestor(v1, v2);

  // Compute gradient contribution of DoFs between v1 and nca
  while (v1 != nca){
    KinEdge* parentEdge = v1->m_parent->findEdge(v1);
    int dof_id = parentEdge->getDOF()->getIndex();
    Math3D::Vector3 deriv = parentEdge->getDOF()->getDerivative(a1->m_position);

    gradient[dof_id] += deriv.dot(diff);
    assert(!std::isnan(gradient[dof_id]));

    v1 = v1->m_parent;
    dofCount++;
  }

  diff = a1->m_position - a2->m_position;
  diff *= (dist - desiredDist) / dist;

  // Compute gradient contribution of DoFs between v2 and nca
  while (v2 != nca){
    KinEdge* parentEdge = v2->m_parent->findEdge(v2);
    int dof_id = parentEdge->getDOF()->getIndex();
    Math3D::Vector3 deriv = parentEdge->getDOF()->getDerivative(a2->m_position);

    gradient[dof_id] += deriv.dot(diff);
    assert( !std::isnan(gradient[dof_id]) );

    v2 = v2->m_parent;
    dofCount++;
  }

  log("debug") << "AtomPairDistanceDirection - DOFs perturbed: " << dofCount << std::endl;
  std::copy(&gradient[0], &gradient[ret->size], ret->data);
}


