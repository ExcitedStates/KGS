//
// Created by Dominik Budday on 03.08.16.
//

#include <cmath>
#include "FixedLink.h"
#include "math3d/primitives.h"
#include "Logger.h"


const double FixedLink::m_maxValue = 0;

Math3D::Vector3 FixedLink::getDerivative(Coordinate& coord) const
{
  return Math3D::Vector3(0,0,0);
}

double FixedLink::getGlobalValue() const
{
  return m_maxValue;
}

double FixedLink::getMaxValue() const
{
  return m_maxValue;
}

double FixedLink::getRandomPerturbation() const
{
  return 0.0;
}

double FixedLink::getMaxPerturbation() const
{
  return m_maxValue;
}

void FixedLink::updateEndVertexTransformation()
{
  m_edge->EndVertex->m_transformation = m_edge->StartVertex->m_transformation;
}