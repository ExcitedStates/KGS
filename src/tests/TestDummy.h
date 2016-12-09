#ifndef TESTDUMMY_H
#define TESTDUMMY_H

#include "TestSuite.h"

#include <string>

using namespace std;

class TestDummy : public TestSuite
{
public:
	bool runTests();
	string name();

private:
	void torsionAccumTest();
	void writeConfiguration();
	void printTorsions();
	void derivativeTest();
    void sampleAndDisplay();
	void sampleHIV1TAR();
};

#endif // TESTDUMMY_H
