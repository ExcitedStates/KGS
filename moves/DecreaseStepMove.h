//
// Created by Rasmus Fonseca on 4/27/16.
//

#ifndef KGS_DECREASESTEPMOVE_H
#define KGS_DECREASESTEPMOVE_H


#include "Move.h"

/**
 * A move which attempts to perform a 'sub-move' a certain amount of steps, while each time
 * decreasing the step-size in the hope that it will eventually be accepted with sufficiently
 * small stepSize.
 */
class DecreaseStepMove: public Move {

 public:
  DecreaseStepMove(Move* m, unsigned int decreaseSteps, double decreaseFactor);

  Configuration* performMove(Configuration* current, gsl_vector* gradient);
 private:
  Move* const m_subMove;
  const unsigned int m_decreaseSteps;
  const double m_decreaseFactor;

};


#endif //KGS_DECREASESTEPMOVE_H
