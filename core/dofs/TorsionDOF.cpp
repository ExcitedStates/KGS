
#include <cmath>
#include "TorsionDOF.h"
#include "math3d/primitives.h"
#include "Logger.h"


const double TorsionDOF::m_maxValue = 0.05;

Math3D::Vector3 TorsionDOF::getDerivative(Coordinate& coord) const
{
  Coordinate& p1 = m_edge->getBond()->Atom1->m_position;
  Coordinate& p2 = m_edge->getBond()->Atom2->m_position;
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

  Coordinate& p1 = m_edge->getBond()->Atom1->m_position;
  Coordinate& p2 = m_edge->getBond()->Atom2->m_position;
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