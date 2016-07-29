
#include <cmath>
#include <cassert>

#include "GlobalRotateDOF.h"
#include "math3d/primitives.h"

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
  //Rotation around origin
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

double GlobalRotateDOF::getMaxValue() const
{
  return 0.01;
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

