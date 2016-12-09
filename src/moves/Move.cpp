#include "Move.h"
#include <math/gsl_helpers.h>
#include <cassert>
#include "Logger.h"

Move::Move():
    m_stepSize(1.0)
{}

Move::~Move(){}

Configuration* Move::move(Configuration* current, gsl_vector* gradient)
{
  return performMove(current, gradient);
}

void Move::setStepSize(double stepSize)
{
  assert(stepSize>0.0);

  m_stepSize = stepSize;
}

double Move::getStepSize()
{
  return m_stepSize;
}
