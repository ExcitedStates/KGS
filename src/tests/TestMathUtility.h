#ifndef TESTMATHUTILITY_H
#define TESTMATHUTILITY_H

#include "TestSuite.h"
#include <string>

using namespace std;

class TestMathUtility : public TestSuite
{
public:
	bool runTests();
	string name();
private:
	bool testFindRotationMatrix();
	bool testTorsionalAngle();
};

#endif // TESTMATHUTILITY_H
