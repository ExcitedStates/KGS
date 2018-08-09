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

#include "Move.h"
#include <math/gsl_helpers.h>
#include <cassert>
#include "Logger.h"

Move::Move():
    m_maxRotation(3.1415/18),//overwrite with getter/setter
    m_scale(false) // by default scaling is disabled
{}

Move::Move(double maxRotation):
    m_maxRotation(maxRotation),//overwrite with getter/setter
    m_scale(true) // by default scaling is enabled, disable via setter if not desired
{}

Move::~Move(){}

Configuration* Move::move(Configuration* current, gsl_vector* gradient)
{
  return performMove(current, gradient);
}

void Move::setMaxRotation(double maxRotation)
{
  assert(maxRotation>0.0);
  m_maxRotation = maxRotation;
}

double Move::getMaxRotation()
{
  return m_maxRotation;
}

void Move::setScalingFlag(bool scale)
{
  m_scale = scale;
}

bool Move::getScalingFlag()
{
  return m_scale;
}