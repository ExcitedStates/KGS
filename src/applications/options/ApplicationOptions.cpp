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

#include "ApplicationOptions.h"

#include "Logger.h"

using namespace std;


ApplicationOptions::ApplicationOptions()
{}

ApplicationOptions::ApplicationOptions(int& argc, char* argv[])
{
  for(int i = 1; i < argc-1; i++) {
    string arg=argv[i];
    if( arg == "--logger" ) {
      string logName = argv[i + 1];
      size_t sep = logName.find(":");
      if(sep != string::npos){
        string logFileName = logName.substr(sep+1);
        logName = logName.substr(0, sep);
        ofstream* os = new ofstream();
        os->open(logFileName.c_str());
        loggerStreams.push_back(os);
        enableLogger(logName, *os);
      }else {
        enableLogger(logName);
      }
    }

    i++;
    //for(int j=i; j<argc-2; j++){
    //  argv[j] = argv[j+2];
    //}
    //argc-=2;
    //i--;
  }
}

/** From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
bool ApplicationOptions::fileExists (const std::string& name)
{
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}
