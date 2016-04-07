
#include "DOF.h"

DOF::DOF(KinEdge* edge):
    m_edge(edge),
    m_value(0)
{ }


void DOF::setValue(double val)
{
  m_value = val;
}

double DOF::getValue() const
{
  return m_value;
}
