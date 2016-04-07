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
#ifndef MATHUTILITY_H
#define MATHUTILITY_H

#include <math.h>
#include <iostream>
#include <math/math.h>

#include "math3d/primitives.h"


#ifndef CTK_PI
#define CTK_PI 3.14159265
#endif

// Trigonometry

double formatRangeRadian (double r); // so that r is in (-Pi,Pi] in radian
Math3D::Matrix3 FindRotationMatrix(const Math3D::Vector3 &axis, double angle); // rotation is according to left hand rule along the axis, NOT right hand rule
Math3D::Vector3 FindRotationArm (Math3D::Vector3 a1, Math3D::Vector3 a2, Math3D::Vector3 p); // a1,a2,p are 3 points in 3D space. Let v=a1-a2. This function returns a Math3D::Vector3 r which is the distance from p to v, that is, r is perpendicular to v, and starts from p.
Math3D::Vector3 ComputeJacobianEntry (Math3D::Vector3 a1, Math3D::Vector3 a2, Math3D::Vector3 p); // compute the jacobian entry (dx,dy,dz) of p when a small rotation is along a1-a2.

/** Compute the torsional angle between the two planes spanned by p1, p2, p3 and p2, p3, p4 respectively. */
double TorsionalAngle (Math3D::Vector3 &p1, Math3D::Vector3 &p2, Math3D::Vector3 &p3, Math3D::Vector3 &p4);
double VectorRotationAngle (Math3D::Vector3 p1, Math3D::Vector3 p2, Math3D::Vector3 axis);
double VectorAngle(Math3D::Vector3 p1, Math3D::Vector3 p2);
double Angle(const Math3D::Vector3& p1, const Math3D::Vector3& p2, const Math3D::Vector3& p3);
double VectorLength (Math3D::Vector3 p1, Math3D::Vector3 p2);
double toDegree(double radian);
double toRadian(double degree);
int Signum(double val);

// Random number generators

/** Generate a random number between 0 and 1 */
double Random01();
/** Generate a random variable between -1 and 1 */
double RandomN1P1();
/** Generate a random normally distributed variable with specified mean and std. deviation. */
double RandomNormalNPiPPi (double mean, double sdv2);
double RandomNormalN1P1 (double mean, double sdv2);
/** Generate a random variable between `-A` and `A` */
double RandomAngleUniform (double A);

#endif
