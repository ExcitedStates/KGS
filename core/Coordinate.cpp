/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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
#include <iostream>
#include <math.h>

#include "Coordinate.h"
#include "Util.h"

using namespace std;

Coordinate::Coordinate () {
	x = 0;
	y = 0;
	z = 0;
}

Coordinate::Coordinate (double X, double Y, double Z) {
	x = X;
	y = Y;
	z = Z;
}

Coordinate::~Coordinate () {
}

string Coordinate::tostring () const {
	string s = Util::d2s(x)+"\t"+Util::d2s(y)+"\t"+Util::d2s(z);
	return s;
}

double Coordinate::distanceTo (Coordinate& other) const {
	return sqrt((x-other.x)*(x-other.x) + (y-other.y)*(y-other.y) + (z-other.z)*(z-other.z));
}

bool Coordinate::isWithinSphere (Coordinate& center, double radius) const {
	return pow(x-center.x,2)+pow(y-center.y,2)+pow(z-center.z,2) <= pow(radius,2);
}

double Coordinate::getAngle (Coordinate& left, Coordinate& right) const {
	double a = distanceTo(left);
	double b = distanceTo(right);
	double c = left.distanceTo(right);
	double angle = acos( ( a*a + b*b - c*c) / (2*a*b) ) * 180 / CTK_PI;
	return angle;
}

Coordinate Coordinate::mid_point (Coordinate& other) const {
	return Coordinate( (x+other.x)/2, (y+other.y)/2, (z+other.z)/2 );
}

Coordinate Coordinate::crossproduct (Coordinate& other) const {
	return Coordinate( y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x );
}

void Coordinate::copyToGslVector (Vector3 v, gsl_vector* gslv) {
	gsl_vector_set(gslv,0,v.x);
	gsl_vector_set(gslv,1,v.y);
	gsl_vector_set(gslv,2,v.z);
}

double Coordinate::getPlanarAngle (Coordinate& c1, Coordinate& c2, Coordinate& c3, Coordinate& c4) {
	Coordinate c21 = Coordinate(c1.x-c2.x,c1.y-c2.y,c1.z-c2.z);
	Coordinate c23 = Coordinate(c3.x-c2.x,c3.y-c2.y,c3.z-c2.z);
	Coordinate c43 = Coordinate(c3.x-c4.x,c3.y-c4.y,c3.z-c4.z);
	Coordinate normal_123 = c21.crossproduct(c23);
	Coordinate normal_234 = c23.crossproduct(c43);
	Coordinate origin(0,0,0);
	double angle = origin.getAngle(normal_123,normal_234);
	return angle;
}
