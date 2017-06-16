//
// Created by Rasmus Fonseca on 6/16/17.
//

#ifndef KGS_APPLICATIONOPTIONS_H
#define KGS_APPLICATIONOPTIONS_H

#include <string>
#include <vector>
#include <fstream>

class ApplicationOptions {

 protected:
  /**
   * Initialize options that are common to all options. Currently only "--logger" arguments are parsed here.
   */
  ApplicationOptions(int& argc, char* argv[]);
  ApplicationOptions();

  bool fileExists (const std::string& name);

 private:
  std::vector<std::ofstream*> loggerStreams;
};


#endif //KGS_APPLICATIONOPTIONS_H
