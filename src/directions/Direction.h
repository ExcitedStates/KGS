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


#ifndef KGS_DIRECTION_H
#define KGS_DIRECTION_H

#include <gsl/gsl_vector.h>

#include "core/Configuration.h"

/**
 * Implementing sub-classes specify different methods of computing gradients.
 */
class Direction {
 public:
  virtual ~Direction() = 0;

  /**
   * Compute a gradient for the configuration `conf` and store the result in `ret`.
   * A target configuration can be specified if the implementing subclass performs directed
   * gradients. If not the target will be ignored.
   *
   * Any values stored in the `ret` vector will be overwritten by this function so there is
   * no need to reset it to zero before use.
   *
   * \pre  { `conf` and `ret` are both non-nullptr }
   * \post { All values in `ret` are guaranteed to be in the the range -π to π. }
   */
  void gradient(Configuration* conf, Configuration* target, gsl_vector* ret);

 protected:
  virtual void computeGradient(Configuration* conf, Configuration* target, gsl_vector* ret) = 0;

};


#endif //KGS_DIRECTION_H
