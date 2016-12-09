#ifndef TESTLOCALREBUILD_H
#define TESTLOCALREBUILD_H

#include <string>

#include "TestSuite.h"

using namespace std;

class TestLocalRebuild : public TestSuite
{
public:
    TestLocalRebuild();
    bool runTests();
    string name();
};

#endif // TESTLOCALREBUILD_H
