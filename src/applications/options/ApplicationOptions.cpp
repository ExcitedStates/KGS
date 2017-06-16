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
