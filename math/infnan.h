/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef MATH_INFNAN_H
#define MATH_INFNAN_H

#include <math.h>

/** @file math/infnan.h
 * @ingroup Math
 * @brief Cross-platform infinity and not-a-number routines.
 *
 * Not necessarily throroughly tested.  Developed partially
 * because cygwin's isnan's go into infinite loops.
 */

namespace Math { 

/** @addtogroup Math */
/*@{*/


#ifndef INFINITY
  #include <limits>
#endif //INFINITY


#ifdef INFINITY
  const static double dInf = INFINITY;
  const static float fInf = INFINITY;
#else
  const static double dInf = std::numeric_limits<double>::infinity();
  const static float fInf = std::numeric_limits<float>::infinity();
#endif // INFINITY


///Returns nonzero if x is not-a-number (NaN)
int IsNaN(double x);
int IsNaN(float x);
///Returns nonzero unless x is infinite or a NaN
int IsFinite(double x);
int IsFinite(float x);
///Returns +1 if x is +inf, -1 if x is -inf, and 0 otherwise
int IsInf(double x);
int IsInf(float x);

/*@}*/
} //namespace Math

#endif
