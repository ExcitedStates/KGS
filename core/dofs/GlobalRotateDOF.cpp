
#include <cmath>
#include <cassert>

#include "GlobalRotateDOF.h"
#include "math3d/primitives.h"

GlobalRotateDOF::GlobalRotateDOF(const KinEdge* edge, int axis):
    DOF(edge),
    m_axis(axis)
{
  assert(axis>=0 && axis<=2);
}

Math3D::Vector3 GlobalRotateDOF::getDerivative(Coordinate& coord) const
{
  Math3D::Vector3 axis(0,0,0);
  axis[m_axis]=1.0;
  return std::move( Math3D::cross( axis, coord ) );
}

double GlobalRotateDOF::getGlobalValue() const
{
  //TODO: Implement
//  std::cerr<<"GlobalRotateDOF::getGlobalValue not implemented"<<std::endl;
  return 0;
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

  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * tr;
}

