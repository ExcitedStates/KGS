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
