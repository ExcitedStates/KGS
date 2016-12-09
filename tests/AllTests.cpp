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

