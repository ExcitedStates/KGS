#ifndef TESTSUGARPUCKER_H
#define TESTSUGARPUCKER_H

#include <string>
#include <iostream>
#include "TestSuite.h"

using namespace std;

class TestSugarPucker : public TestSuite
{
public:
    string name(){ return "Sugar puckering"; }
    bool runTests();
private:
    bool testProtein();
    bool testRNA();
    bool testTransformationPropagation();
};

#endif // TESTSUGARPUCKER_H
