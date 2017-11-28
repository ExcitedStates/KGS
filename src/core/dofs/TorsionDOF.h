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
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_TORSIONDOF_H
#define KGS_TORSIONDOF_H


#include "DOF.h"

class TorsionDOF: public DOF {
 public:
  TorsionDOF(const KinEdge* edge): DOF(edge){}
  ~TorsionDOF(){};

  Math3D::Vector3 getDerivative(Coordinate& coord) const override;

  double getGlobalValue() const override;

  double getRandomPerturbation() const override;
  double getMaxPerturbation() const override;

 protected:

  void updateEndVertexTransformation() override;

 private:
  static const double m_maxValue;
};


#endif //KGS_TORSIONDOF_H
