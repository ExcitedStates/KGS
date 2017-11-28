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


#ifndef KGS_LSNRELATIVEDIRECTION_H
#define KGS_LSNRELATIVEDIRECTION_H

#include <gsl/gsl_vector.h>

#include "Selection.h"
#include "Direction.h"
#include "core/Configuration.h"

class LSNrelativeDirection: public Direction  {
 public:
  LSNrelativeDirection(Selection& atomsMoving,
                       const std::vector< std::tuple<Atom*, Atom*, double> >& goal_distances);

 protected:
  void computeGradient(Configuration* conf, Configuration* conf2, gsl_vector* ret) override;

 private:

  void fillmatrices(Configuration* current_q, gsl_matrix* targetJacobian, gsl_matrix* targetDirection);
  // void clashFreeGradient(gsl_vector* gradient, gsl_vector* admissible_gradient, Molecule* protein);
  gsl_matrix* determineBestMove(gsl_matrix* N, gsl_matrix* targetJacobian, gsl_matrix* TargetPosition);
  Selection& m_atomsMovingSelection;
  std::vector< std::tuple<Atom*, Atom*, double> > goal_distances;
};



#endif //KGS_LSNRELATIVEDIRECTION_H
