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

//
// Created by Rasmus Fonseca on 4/27/16.
//

#include <math/gsl_helpers.h>
#include "DecreaseStepMove.h"
#include "core/Molecule.h"

DecreaseStepMove::DecreaseStepMove(Move* m, unsigned int decreaseSteps, double decreaseFactor):
    Move(),
    m_subMove(m),
    m_decreaseSteps(decreaseSteps),
    m_decreaseFactor(decreaseFactor)
{
}

Configuration* DecreaseStepMove::performMove(Configuration* current, gsl_vector* gradient) {
  double subStepSize = m_stepSize;

  m_subMove->setStepSize(subStepSize);
  Configuration* new_q;

  //As sub-moves sometime modify the gradient we need a copy
  gsl_vector* gradient_copy = gsl_vector_copy(gradient);


  //If resulting structure is in collision try scaling down the gradient
  for (int step = 0; step <= m_decreaseSteps; step++) {
    new_q = m_subMove->move(current, gradient_copy);

    if (new_q->updatedMolecule()->inCollision()) {
      gsl_vector_memcpy(gradient_copy, gradient);
      delete new_q;

      subStepSize *= m_decreaseFactor;
      m_subMove->setStepSize(subStepSize);

    } else {
      gsl_vector_free(gradient_copy);
      return new_q;
    }
  }

  gsl_vector_free(gradient_copy);
  return new_q;
}
