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

#include "TestMathUtility.h"
#include <iomanip>
#include "../Logger.h"
#include "math/MathUtility.h"
#include "../utils/math3d/primitives.h"

bool TestMathUtility::runTests(){
    if(testFindRotationMatrix()) log("test")<<left<<setw(60)<<"TestMathUtility::testFindRotationMatrix:"<<"passed"<<endl;
    else { log("test")<<left<<setw(60)<<"TestMathUtility::testFindRotationMatrix:"<<"failed"<<endl;return false;}
    if(testTorsionalAngle()) log("test")<<left<<setw(60)<<"TestMathUtility::testTorsionalAngle:"<<"passed"<<endl;
    else { log("test")<<left<<setw(60)<<"TestMathUtility::testTorsionalAngle:"<<"failed"<<endl;return false;}
    return true;

}

bool TestMathUtility::testFindRotationMatrix(){
	Vector3 axis(1,1,1);
	double angle = 3.141592*2.0/3;
	Matrix3 m = MathUtility::FindRotationMatrix(axis, angle);
	Vector3 p(3,0,0);
	p = m*p;
	cout<<p<<endl;
	if( fabs(p.x-0)>0.0001 || fabs(p.y-0)>0.0001 || fabs(p.z-3)>0.0001 ){
		log("test")<<"TestMathUtility::testFindRotationMatrix(): Expected (0,0,3) but got "<<p<<endl;
		return false;
	}


	return true;
}
bool TestMathUtility::testTorsionalAngle(){
	Vector3 p1(0,-1,0);
	Vector3 p0(0,0,0);
	Vector3 axis(1,1,1);
	double angle = 3.141592*2.0/3;
	Matrix3 m = MathUtility::FindRotationMatrix(axis, angle);
	Vector3 p(3,0,0);
	double before = TorsionalAngle(p1,p0,axis,p);
	//cout<<"TorsionalAngle("<<p1<<","<<p0<<","<<axis<<","<<p<<") = "<<before<<endl;
	p = m*p;
	double after = TorsionalAngle(p1,p0,axis,p);
	//cout<<"TorsionalAngle("<<p1<<","<<p0<<","<<axis<<","<<p<<") = "<<after<<endl;
	if(fabs((after-before)+angle)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 1: Expected "<<(-angle)<<" but got "<<(after-before)<<endl;
		return false;
	}

	Vector3 p2, p3;

	p0.x = 0; p0.y = 1; p0.z = 0;
	p1.x = 1; p1.y = 1; p1.z = 1;
	p2.x = 1; p2.y = 2; p2.z = 1;
	p3.x = 3; p3.y = 2; p3.z = 3;
	double torsion = TorsionalAngle(p0,p1,p2,p3);
	double expected = CTK_PI;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 2: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}
	p2.x = 1.1; p2.y = 3; p2.z = 1.1;
	torsion = TorsionalAngle(p0,p1,p2,p3);
	expected = CTK_PI;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 3: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}

	p1.x = 1; p1.y = 1; p1.z = 0;
	p2.x = 1; p2.y = 2; p2.z = 0;
	p3.x = 1; p3.y = 2; p3.z = 1;
	torsion = TorsionalAngle(p0,p1,p2,p3);
	expected = CTK_PI/2;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 4: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}

	p3.x = 1+1/sqrt(2); p3.y = 2; p3.z = 1/sqrt(2);
	torsion = TorsionalAngle(p0,p1,p2,p3);
	expected = 3*CTK_PI/4;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 5: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}

	p3.x = 1+1/sqrt(2); p3.y = 2; p3.z = -1/sqrt(2);
	torsion = TorsionalAngle(p0,p1,p2,p3);
	expected = -3*CTK_PI/4;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 6: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}

	p0.x = -1.989; p0.y = -9.046; p0.z = 6.143;
	p1.x = -0.574; p1.y = -8.494; p1.z = 5.280;
	p2.x = -0.719; p2.y = -8.681; p2.z = 3.704;
	p3.x = -0.490; p3.y = -9.948; p3.z = 3.094;
	torsion = TorsionalAngle(p0,p1,p2,p3);
	expected = -1.40247;
	if(fabs(torsion-expected)>0.0001){
		log("test")<<"TestMathUtility::testTorsionalAngle() 7: Expected "<<expected<<" but got "<<torsion<<endl;
		return false;
	}
	return true;
}

string TestMathUtility::name(){
	return "MathUtility";
}
