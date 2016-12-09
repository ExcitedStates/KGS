#include <cassert>
#include "DOF.h"

DOF::DOF(const KinEdge* edge):
    m_edge(edge),
    m_value(0),
    m_index(-1),
    m_cycleIndex(-1)
{ }


void DOF::setValue(double val)
{
  m_value = val;
}

double DOF::getValue() const
{
  return m_value;
}

unsigned int DOF::getIndex() const
{
  return (unsigned int)m_index;
}

unsigned int DOF::getCycleIndex() const
{
  return (unsigned int)m_cycleIndex;
}

void DOF::setIndex(unsigned int idx)
{
  m_index = idx;
}

void DOF::setCycleIndex(unsigned int idx)
{
  m_cycleIndex = idx;
}
