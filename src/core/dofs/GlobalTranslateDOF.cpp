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


#include <cmath>
#include <cassert>

#include "GlobalTranslateDOF.h"
#include "math3d/primitives.h"
#include <iostream>

const double GlobalTranslateDOF::m_maxValue = 1.0;

GlobalTranslateDOF::GlobalTranslateDOF(const KinEdge* edge, int axis):
    DOF(edge),
    m_axis(axis)
{
  assert(axis>=0 && axis<=2);

}

Math3D::Vector3 GlobalTranslateDOF::getDerivative(Coordinate& coord) const
{
  Math3D::Vector3 ret(0,0,0);
  ret[m_axis]=1.0;
  return ret;
}

double GlobalTranslateDOF::getGlobalValue() const
{
  //TODO: Implement
//  std::cerr<<"GlobalTranslateDOF::getGlobalValue not implemented"<<std::endl;
  return 0;
}

double GlobalTranslateDOF::getRandomPerturbation() const
{
  return m_maxValue*RandomN1P1();
}

double GlobalTranslateDOF::getMaxPerturbation() const
{
  return m_maxValue;
}

void GlobalTranslateDOF::updateEndVertexTransformation()
{
  //if( std::fabs(m_value)<0.000001 ) {
  //  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation;
  //  return;
  //}

  Math3D::RigidTransform tr;
  tr.setIdentity();
  tr.t.data[m_axis]=m_value;

  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * tr;
}
