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




#include <gsl/gsl_vector_double.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "math/gsl_helpers.h"
#include "BlendedDirection.h"
#include <iostream>

BlendedDirection::BlendedDirection():
    m_totalWeight(0)
{ }

void BlendedDirection::addDirection(Direction* dir, double weight)
{
  m_directions.push_back(dir);
  m_weights.push_back(weight);
  m_totalWeight+=weight;
}

void BlendedDirection::changeWeight(int dirCount, double weight)
{
      m_totalWeight -= m_weights[dirCount];
      m_weights[dirCount] = weight;
      m_totalWeight += weight;
}

void BlendedDirection::computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret)
{
  assert(!m_directions.empty());

  m_directions[0]->gradient(conf,target,ret);
  double magnitude = gsl_vector_length(ret);

  gsl_vector* tmp = gsl_vector_calloc(ret->size);
  for(size_t i=1;i<m_directions.size();i++){
    m_directions[i]->gradient(conf,target,tmp);
    gsl_vector_scale(tmp,magnitude);

    for (int entry=0; entry < ret->size; ++entry){
      double val = m_weights[i] * gsl_vector_get(tmp,entry) + m_weights[0] * gsl_vector_get(ret,entry);
      //TODO, this doesn't work or more than 2 directions.
//      gsl_vector_set(ret,i,formatRangeRadian(val));
      gsl_vector_set(ret,i,val);
    }
  }
}
