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

#include <vector>
#include <math/gsl_helpers.h>
#include <cmath>

#include "RandomDirection.h"
#include "SamplingOptions.h"
#include "core/RigidbodyGraph.h"
#include "core/Molecule.h"
#include "math.h" //dPi
#include "math/MathUtility.h"

RandomDirection::RandomDirection(double maxRotation):
    m_maxRotation(maxRotation)
{ }


void RandomDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  //TODO: It shouldn't be necessary to go through edges of the tree ...
  //TODO: Just set all entries in ret to random values (?)
  Molecule * protein = conf->getProtein();
  double absMax = 0.0;
  for (auto const& edge: protein->m_spanning_tree->Edges) {
    int dofId = edge->DOF_id;
    double newVal = RandomAngleUniform(dPi);
    absMax = std::max(absMax, std::fabs(newVal));

    //ToDo: adapt to sampling options selectionMoving etc
    const std::vector<int> &res_network = SamplingOptions::getOptions()->residueNetwork;
    if (res_network.empty()) {
      gsl_vector_set(ret, dofId, newVal);
    } else {
      int resId = edge->getBond()->Atom1->getResidue()->getId();
      if (std::find(res_network.begin(), res_network.end(), resId) != res_network.end()) {
        gsl_vector_set(ret, dofId, newVal);
      }
    }
  }

  if ( absMax > m_maxRotation ){
    gsl_vector_scale(ret, m_maxRotation/absMax);
  }

  gsl_vector_scale(ret, SamplingOptions::getOptions()->stepSize);
}
