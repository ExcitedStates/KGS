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



#include <vector>
#include <math/gsl_helpers.h>
#include <cmath>

#include "RandomDirection.h"

RandomDirection::RandomDirection(Selection& selectionMoving, double maxRotation):
    m_maxRotation(maxRotation),
    m_selectionMoving(selectionMoving)
{
}


void RandomDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  //TODO: It shouldn't be necessary to go through m_edges of the tree ...
  //TODO: Just set all entries in ret to random values (?)
  Molecule * protein = conf->getMolecule();
  double absMax = 0.0;
  for (auto const& edge: protein->m_spanningTree->m_edges) {
    int dofId = edge->getDOF()->getIndex();
    double newVal = 0.0;

    if( edge->getBond()!=nullptr && m_selectionMoving.inSelection(edge->getBond())) {
      newVal = edge->getDOF()->getRandomPerturbation();
    }

    absMax = std::max(absMax, std::fabs(newVal));

    gsl_vector_set(ret, dofId, newVal);
  }

  if ( absMax > m_maxRotation ){
    gsl_vector_scale(ret, m_maxRotation/absMax);
  }
}
