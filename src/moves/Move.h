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

  double getMaxRotation();
  void setMaxRotation(double maxRotation);

  virtual void setScalingFlag(bool scale); //may be overridden in composite moves etc.
  bool getScalingFlag();

 protected:
  Move();
  Move(double maxRotation);

  virtual Configuration* performMove(Configuration* current, gsl_vector* gradient) = 0;

  int m_movesRejected;
  int m_movesAccepted;
  double m_maxRotation;
  bool m_scale;
};

#endif

