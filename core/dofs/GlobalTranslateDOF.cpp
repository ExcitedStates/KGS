
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