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
