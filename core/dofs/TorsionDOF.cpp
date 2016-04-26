
#include <cmath>
#include "TorsionDOF.h"
#include "math3d/primitives.h"

Math3D::Vector3 TorsionDOF::getDerivative(Coordinate& coord) const
{
  Coordinate& p1 = m_edge->getBond()->Atom1->m_Position;
  Coordinate& p2 = m_edge->getBond()->Atom2->m_Position;
  Math3D::Vector3 axis = p2-p1;
  axis.inplaceNormalize();
  Math3D::Vector3 arm = coord-p1;
  return std::move( Math3D::cross( axis, arm ) );
}

double TorsionDOF::getGlobalValue() const
{
  //TODO: Implement
  return 0;
}

void TorsionDOF::updateEndVertexTransformation()
{
  //if( std::fabs(m_value)<0.000001 ) {
  //  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation;
  //  return;
  //}

  Coordinate& p1 = m_edge->getBond()->Atom1->m_Position;
  Coordinate& p2 = m_edge->getBond()->Atom2->m_Position;
  Math3D::Vector3 axis = p2-p1;
  axis.inplaceNormalize();

  Math3D::RigidTransform m1, m2, m3;
  m1.setIdentity();
  m2.setIdentity();
  m3.setIdentity();

  //TODO: Theres a closed-form expression that would save some matrix multiplications and make this faster
  m1.setTranslate(p1);
  m2.setRotate(FindRotationMatrix(axis, -m_value)); //FindRotationMatrix returns left-handed rotation
  m3.setTranslate(-1.0*p1);

  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation * m1*m2*m3;
}
