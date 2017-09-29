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
#include "Logger.h"

using namespace std;

Logger::Logger() {}

void Logger::enableLogger(const string& name, ostream& stream){
  activeLoggers[name] = &stream;
}

void Logger::disableLogger(const string& name){
  activeLoggers.erase(name);
}


ostream& Logger::log(const string& name){
  if( !loggerEnabled(name) )
    return drain;

  return *(activeLoggers[name]);
}

Logger* Logger::getInstance(){
  if(Logger::instance==nullptr) Logger::instance = new Logger();
  return Logger::instance;
}

bool Logger::loggerEnabled(const string& name){
  return activeLoggers.count(name)!=0;
}


Logger* Logger::instance = nullptr;


ostream& log(){ return Logger::getInstance()->log("default"); }

ostream& log(const string& loggerName){ return Logger::getInstance()->log(loggerName); }

void enableLogger(const string& loggerName) { Logger::getInstance()->enableLogger(loggerName, cout); }

void enableLogger(const string& loggerName, ostream& stream) { Logger::getInstance()->enableLogger(loggerName, stream); }

void disableLogger(const string& loggerName) { Logger::getInstance()->disableLogger(loggerName); }

bool loggerEnabled(const string& loggerName) { return Logger::getInstance()->loggerEnabled(loggerName); }


