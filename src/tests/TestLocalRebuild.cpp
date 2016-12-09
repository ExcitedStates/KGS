#include "TestLocalRebuild.h"

#include <iostream>

#include "../Logger.h"

TestLocalRebuild::TestLocalRebuild()
{
}

bool TestLocalRebuild::runTests(){
    log("test")<<"TestLocalRebuild: No tests yet. FAILED."<<endl;
    return false;
}

string TestLocalRebuild::name(){
    return "LocalRebuild";
}
