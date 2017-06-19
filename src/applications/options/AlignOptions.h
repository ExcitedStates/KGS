//
// Created by Dominik Budday on 16.06.17.
//

#ifndef KGS_ALIGNOPTIONS_H
#define KGS_ALIGNOPTIONS_H

#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"

class Molecule;

class Atom;

class AlignOptions {
 public:

  /** The working directory */
  std::string workingDirectory;
  /** File-path for the initial structure. */
  std::string initialStructureFile;
  /** File-path for the target structure. */
  std::string targetStructureFile;

  /** Specifies the residues of the molecule that will undergo RMSD alignment during sampling. */
  std::string alignSelection;

  void print();

 private:
  AlignOptions();

  AlignOptions(int argc, char *argv[]);

  void printUsage(char *pname);

  inline bool fileExists(const std::string &name);

  void initializeVariables();

  //Singleton pattern
  static AlignOptions *instance;
 public:
  static AlignOptions *getOptions();

  static AlignOptions *createOptions(int argc, char *argv[]);

  static AlignOptions *createOptions();

};

#endif //KGS_ALIGNOPTIONS_H
