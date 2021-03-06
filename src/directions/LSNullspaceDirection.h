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


#ifndef KGS_LSnullptrSPACEDIRECTION_H
#define KGS_LSnullptrSPACEDIRECTION_H

#include <gsl/gsl_vector.h>

#include "Selection.h"
#include "Direction.h"
#include "core/Configuration.h"


class LSNullspaceDirection: public Direction {
public:
  LSNullspaceDirection(Selection& atomsMoving);

private:
  void computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret);

  void fillmatrices(Configuration* current_q, Configuration* pTarget, gsl_matrix* targetJacobian, gsl_matrix* targetPosition);
  void clashFreeGradient(gsl_vector* gradient, gsl_vector* admissible_gradient, Molecule* protein);
  gsl_matrix* determineBestMove(gsl_matrix* N, gsl_matrix* targetJacobian, gsl_matrix* TargetPosition);


//  static gsl_matrix* m_TargetJacobian;
//  static gsl_matrix* m_TargetPosition;

  Selection& m_atomsMovingSelection;
};


#endif //KGS_LSnullptrSPACEDIRECTION_H
