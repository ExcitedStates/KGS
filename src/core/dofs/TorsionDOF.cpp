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
#include "TorsionDOF.h"
#include "math3d/primitives.h"
#include "Logger.h"


const double TorsionDOF::m_maxValue = 0.05;

Math3D::Vector3 TorsionDOF::getDerivative(Coordinate& coord) const
{
  Coordinate& p1 = m_edge->getBond()->m_atom1->m_position;
  Coordinate& p2 = m_edge->getBond()->m_atom2->m_position;
  Math3D::Vector3 axis = p2-p1;
  axis.inplaceNormalize();
  Math3D::Vector3 arm = coord-p1;
  return std::move( Math3D::cross( axis, arm ) );
}

double TorsionDOF::getGlobalValue() const
{
  return m_edge->getBond()->getTorsion();
}

double TorsionDOF::getRandomPerturbation() const
{
  return RandomAngleUniform(m_maxValue);
}

double TorsionDOF::getMaxPerturbation() const
{
  return m_maxValue;
}

void TorsionDOF::updateEndVertexTransformation()
{
  //if( std::fabs(m_value)<0.000001 ) {
  //  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation;
  //  return;
  //}

  Coordinate& p1 = m_edge->getBond()->m_atom1->m_position;
  Coordinate& p2 = m_edge->getBond()->m_atom2->m_position;
  Math3D::Vector3 axis = p2-p1;
  axis.inplaceNormalize();

  /// Previous implementation
//  Math3D::RigidTransform m1, m2, m3;
//  m1.setIdentity();
//  m2.setIdentity();
//  m3.setIdentity();
//
//  //Done: Theres a closed-form expression that would save some matrix multiplications and make this faster
//  m1.setTranslate(p1);
//  m2.setRotate(FindRotationMatrix(axis, -m_value)); //FindRotationMatrix returns left-handed rotation
//  m3.setTranslate(-1.0*p1);
//
//  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * m1*m2*m3;

  ///Closed-form expression
  Math3D::RigidTransform m;
  m.setIdentity();
  m.setRotation(FindRotationMatrix(axis, -m_value));
  m.setTranslation(p1 - m.R * p1 );

  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * m;
}
