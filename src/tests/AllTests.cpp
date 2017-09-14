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

#include <iomanip>

#include "AllTests.h"
#include "TestLocalRebuild.h"
#include "TestSugarPucker.h"
#include "TestMathUtility.h"
#include "TestDummy.h"
#include "../Logger.h"
#include <string>

using namespace std;

AllTests::AllTests(){
#ifdef RUNDUMMYTEST
    allTests.push_back(new TestDummy());
#endif
    allTests.push_back(new TestLocalRebuild());
    allTests.push_back(new TestSugarPucker());
    allTests.push_back(new TestMathUtility());
}

string AllTests::name(){ return "All tests"; }

bool AllTests::runTests(){
    enableLogger("test");
    bool allPassed = true;
    for(vector<TestSuite*>::iterator tIt = allTests.begin(); tIt!=allTests.end(); tIt++){
        TestSuite* ts = *tIt;
        bool passed = ts->runTests();
        cout.flush();
        cerr.flush();
        log("test")<<std::left<<std::setw(60)<<ts->name()<<" : "<<(passed?"PASSED":"FAILED")<<endl;
        allPassed = passed && allPassed;
    }
    log("test")<<std::left<<std::setw(60)<<name()<<" : "<<(allPassed?"PASSED":"FAILED")<<endl;
    return allPassed;
}

