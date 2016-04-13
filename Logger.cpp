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


