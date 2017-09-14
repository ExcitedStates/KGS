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

#include "GlobalRotateDOF.h"
#include "math3d/primitives.h"

const double GlobalRotateDOF::m_maxValue = 0.05;

GlobalRotateDOF::GlobalRotateDOF(const KinEdge* edge, int axis):
    DOF(edge),
    m_axis(axis)
{
  assert(axis>=0 && axis<=2);

  KinVertex* v = edge->EndVertex;
  while(v->m_rigidbody==nullptr){
    v = v->m_edges[0]->EndVertex;
  }
  m_firstAtom = v->m_rigidbody->Atoms[0];
}

Math3D::Vector3 GlobalRotateDOF::getDerivative(Coordinate& coord) const
{
  Math3D::Vector3 axis(0,0,0);
  axis[m_axis]=1.0;
//  Rotation around origin
//  return Math3D::cross( axis, coord ) ;

  //Rotation around first atom
  Math3D::Vector3 arm = coord - (m_firstAtom->m_position);
  return Math3D::cross( axis, arm ) ;
}

double GlobalRotateDOF::getGlobalValue() const
{
  //TODO: Implement
//  std::cerr<<"GlobalRotateDOF::getGlobalValue not implemented"<<std::endl;
  return 0;
}


double GlobalRotateDOF::getRandomPerturbation() const
{
  return m_maxValue*RandomAngleUniform(m_maxValue);
}

double GlobalRotateDOF::getMaxPerturbation() const
{
  return m_maxValue;
}

void GlobalRotateDOF::updateEndVertexTransformation()
{
  //if( std::fabs(m_value)<0.000001 ) {
  //  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation;
  //  return;
  //}

  Math3D::RigidTransform tr;
  tr.setIdentity();
  switch(m_axis) {
    case 0: tr.R.setRotateX(m_value); break;
    case 1: tr.R.setRotateY(m_value); break;
    case 2: tr.R.setRotateZ(m_value); break;
  }


  //Rotation around origin
//  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * tr;

  //Rotation around first atom
  Math3D::RigidTransform m1, m3;
  m1.setIdentity();
  m3.setIdentity();
  m1.setTranslate(     m_firstAtom->m_position);
  m3.setTranslate(-1.0*m_firstAtom->m_position);
  m_edge->EndVertex->m_transformation =
      m_edge->StartVertex->m_transformation * m1 * tr * m3;

}