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
