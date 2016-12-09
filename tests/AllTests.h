#ifndef ALLTESTS_H
#define ALLTESTS_H

#include <iostream>
#include <string>
#include <vector>

#include "TestSuite.h"

using namespace std;

class AllTests : public TestSuite
{
public:
    AllTests();
    string name();
    bool runTests();
private:
    vector<TestSuite*> allTests;
};

#endif // ALLTESTS_H
