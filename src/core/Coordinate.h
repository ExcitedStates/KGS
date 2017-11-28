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

#ifndef COORDINATE_H
#define COORDINATE_H

#include <string>

#include "math3d/primitives.h"
#include <gsl/gsl_vector.h>

class Coordinate : public Math3D::Vector3 {
  public:
	Coordinate();
	Coordinate(double x,double y,double z);
	~Coordinate();

	std::string tostring() const;
	double distanceTo (Coordinate& other) const; // Euclidean distance between self and other
	bool isWithinSphere (Coordinate& center, double radius) const;
	double getAngle (Coordinate& left, Coordinate& right) const; // Angle in degrees between left-this-right
	Coordinate mid_point (Coordinate& other) const;
	Coordinate crossproduct (Coordinate& other) const;

	static void copyToGslVector (Math3D::Vector3 v, gsl_vector* gslv);
	static double getPlanarAngle (Coordinate& c1, Coordinate& c2, Coordinate& c3, Coordinate& c4); // Angle between plane c123 and plane c234
};

#endif
