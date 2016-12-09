#ifndef MOVE_H_
#define MOVE_H_

#include "core/Configuration.h"

class Move{
 public:

  virtual ~Move() = 0;

  /**
   * Perform a perturbation of `current` in the direction of the specified gradient.
   * Subclasses of `Move` are permitted to disregard or modify the gradient to satisfy
   * constraints for example, but they will aim to update all DOFs according to the
   * gradient whenever possible.
   *
   * A configuration is always returned, but it is not guaranteed to be clash-free.
   */
  Configuration* move(Configuration* current, gsl_vector* gradient);

  double getStepSize();

  void setStepSize(double stepSize);


 protected:
  Move();

  virtual Configuration* performMove(Configuration* current, gsl_vector* gradient) = 0;

  int m_movesRejected;
  int m_movesAccepted;
  double m_stepSize;
};

#endif

