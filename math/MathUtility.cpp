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
#include "MathUtility.h"

using namespace Math3D;

double formatRangeRadian (double r) {
  while (r<=-dPi) r += 2*dPi;
  while (r>dPi)   r -= 2*dPi;
  return r;
}

Matrix3 FindRotationMatrix(const Vector3 &axis, double angle){
  Matrix3 rotMat;
  Vector3 scaledAxis(axis);
  double c = cos(angle), s = sin(angle), omc = 1 - c, len = scaledAxis.length();

  if (len != 1) scaledAxis /= len;

  double x = scaledAxis[0];
  double y = scaledAxis[1];
  double z = scaledAxis[2];
  double xs = x * s;
  double ys = y * s;
  double zs = z * s;
  double xyomc = x * y * omc;
  double xzomc = x * z * omc;
  double yzomc = y * z * omc;

  rotMat(0,0) = x*x*omc + c;
  rotMat(0,1) = xyomc + zs;
  rotMat(0,2) = xzomc - ys;

  rotMat(1,0) = xyomc - zs;
  rotMat(1,1) = y*y*omc + c;
  rotMat(1,2) = yzomc + xs;

  rotMat(2,0) = xzomc + ys;
  rotMat(2,1) = yzomc - xs;
  rotMat(2,2) = z*z*omc + c;

  return rotMat;
}

Vector3 FindRotationArm (Vector3 a1, Vector3 a2, Vector3 p) {
	Vector3 pa2 = p - a2;
	Vector3 v = a1 - a2;
	Vector3 anchor = a2 + (pa2.dot(v)/(v.norm()*v.norm()))*v;
	return p - anchor;
}

Vector3 ComputeJacobianEntry (Vector3 a1, Vector3 a2, Vector3 p) {
//	Vector3 rotation_arm = FindRotationArm(a1,a2,p); // no need to compute the perpendicular arm
	Vector3 rotation_arm = p - a2;
	Vector3 axis = a2 - a1;
	axis.setNormalized(axis);
	return cross(axis,rotation_arm);
}


//Generate a random double between 0 and 1 */
double Random01()
{
	return (1.0 * rand()) / RAND_MAX;
}

//Generate a random variable between -1 and 1
double RandomN1P1()
{
	return (2.0 * rand()) / RAND_MAX - 1;
}

/* Generate a random angle between -Pi and Pi according to the normal distribution with mean and sdv2 */
double RandomNormalNPiPPi (double mean, double sdv2) {
	double random, p_random, p_control;
	while (true) {
		random = RandomN1P1()*dPi;
		p_random = exp(-0.5*pow((random-mean),2)/sdv2)/(sqrt(2*dPi*sdv2));
		p_control = Random01();
		if (p_random>=p_control) {
			break;
		}
	}
	return random;
}

// generate a random number between -1 and 1 according to Normal distribution with mean and sdv2
double RandomNormalN1P1 (double mean, double sdv2) {
  double random, p_random, p_control;
  while (true) {
    random = RandomN1P1();
    p_random = exp(-0.5*pow((random-mean),2)/sdv2)/(sqrt(2*dPi*sdv2));
    p_control = Random01();
    if (p_random>=p_control) {
      break;
    }
  }
  return random;
}

// generate a random angle between -A and A according to uniform distribution
double RandomAngleUniform (double A) {
	return RandomN1P1()*A;
}

// compute the torsional angle p1-p2-p3-p4 in radians
double TorsionalAngle (Vector3 &p1, Vector3 &p2, Vector3 &p3, Vector3 &p4) {
	Vector3 p21 = p1-p2;
	Vector3 p23 = p3-p2;
	Vector3 p43 = p3-p4;
	Vector3 normal_213 = cross(p21,p23);
	Vector3 normal_432 = cross(p23,p43);
	normalize(normal_213);
	normalize(normal_432);
	normalize(p23);
	double angle = VectorRotationAngle(normal_213,normal_432,p23);
	return angle;
}

double VectorRotationAngle (Vector3 p1, Vector3 p2, Vector3 axis) { 
// The function computes the angle to rotate p1 to p2 around axis clockwisely.
// p1 and p2 are perpendicular to the axis.
// All 3 input vectors are unit vectors
	double cos = dot(p1,p2);
	if (cos<-1.0) {
		cos = -1.0;
	}
	else if (cos>1.0) {
		cos = 1.0;
	}
	double angle = acos(cos);
	double s = axis.x*(p1.z*p2.y-p1.y*p2.z)+axis.y*(p1.x*p2.z-p1.z*p2.x)+axis.z*(p1.y*p2.x-p1.x*p2.y);
	if (s>0.0) {
		angle *= -1;
	}
	return angle;
}

double VectorAngle (Vector3 p1, Vector3 p2) {
  normalize(p1);
  normalize(p2);
  double cos = dot(p1,p2);
  return acos(cos);
}

double Angle(const Math3D::Vector3& p1, const Math3D::Vector3& p2, const Math3D::Vector3& p3)
{
	return VectorAngle(p1-p2, p3-p2);
}

double VectorLength (Vector3 p1, Vector3 p2) {
	Vector3 v = p1-p2;
	return sqrt(dot(v,v));

}

double toDegree(double radian) {
	return radian*dRtoDConst;
}

double toRadian(double degree) {
	return degree*dDtoRConst;
}

int Signum(double val){
	return (val>0) - (val<0);
}
