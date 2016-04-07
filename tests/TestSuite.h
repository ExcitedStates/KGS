#ifndef TESTSUITE_H
#define TESTSUITE_H

#include <string>

using namespace std;

class TestSuite{
public:
    virtual bool runTests() = 0;
    virtual string name() = 0;
};

#endif // TESTSUITE_H
