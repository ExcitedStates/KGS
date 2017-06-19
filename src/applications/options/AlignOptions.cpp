//
// Created by Dominik Budday on 16.06.17.
//

#include "AlignOptions.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdexcept>


#include "Util.h"
#include "Logger.h"
#include "Selection.h"
#include "core/Molecule.h"
#include "core/Atom.h"
#include "core/Chain.h"


using namespace std;

AlignOptions::AlignOptions(){
  initializeVariables();
}

AlignOptions::AlignOptions(int argc, char* argv[])
{
  initializeVariables();

  if(argc<2){
    enableLogger("so");
    printUsage(argv[0]);
    exit(-1);
  }

  int i;
  for(i=1;i<argc;i++){
    string arg = argv[i];
    if(arg=="--initial"){                       initialStructureFile = argv[++i];                   continue; }
    if(arg=="--target"){                        targetStructureFile = argv[++i];                    continue; }
    if(arg=="--workingDirectory"){              workingDirectory = argv[++i];                       continue; }
    if(arg=="--align" or arg=="--alignSelection"){  alignSelection = argv[++i];                     continue; }

    if(arg.at(0)=='-'){
      cerr<<"Unknown option: "<<arg<<endl<<endl;
      enableLogger("so");
      printUsage(argv[0]);
      exit(-1);
    }
  }

  //Check initial structure
  if(initialStructureFile.empty()) {
    enableLogger("so");
    cerr<<"Initial structure file must be supplied"<<endl<<endl;
    exit(-1);
  }
  //Check target structure
  if(targetStructureFile.empty()){
    log("so")<<"Target structure file must be supplied."<<endl<<endl;
    exit(-1);
  }

  //Set workingDirectory and moleculeName using the initialStructureFile.
  char* tmp = realpath(initialStructureFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<initialStructureFile<<" is not a valid PDB-file"<<endl; exit(-1); }
  string pdb_file(tmp);
  int nameSplit = pdb_file.find_last_of("/\\");
  if(workingDirectory.empty()) {
    workingDirectory = pdb_file.substr(0, nameSplit + 1);
    log("so")<<"Changing working directory to "<<workingDirectory<<" using initial"<<endl;
  }

  //Check target structure
  if(targetStructureFile.empty()){
    log("so")<<"No target structure file supplied, random sampling."<<endl;
  }

  free(tmp);
}


void AlignOptions::initializeVariables(){
  initialStructureFile      = "";
  targetStructureFile       = "";
  workingDirectory          = "";
  alignSelection            = "all";
}

void AlignOptions::print(){
  log("so")<<"Sampling options:"<<std::endl;
  log("so")<<"\t--initial "<<initialStructureFile<<endl;
  log("so")<<"\t--target "<<targetStructureFile<<endl;
  log("so")<<"\t--workingDirectory "<<workingDirectory<<endl;
  log("so")<<"\t--align "<<alignSelection<<endl;
}

void AlignOptions::printUsage(char* pname){
  log("so")<<"Usage: "<<pname<<" --initial <pdb-file> --target <pdb-file> [option list]"<<endl;
  log("so")<<"Computes aligned structure and rmsd given specified selections."<<endl;
  log("so")<<endl;
  log("so")<<"Options:"<<endl;

  log("so")<<"  --initial <pdb-file> \t: Specifies the initial structure."<<endl;
  log("so")<<"  --target <pdb-file> \t: Specifies the target structure, optional."<<endl;
  log("so")<<"  --workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;
  log("so")<<"  --align < A pymol-like pattern that indicates which subset of atoms are used during alignment (if specified). Default is 'all'."<<endl;
}



/** From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
inline bool AlignOptions::fileExists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

AlignOptions* AlignOptions::instance = nullptr;

AlignOptions* AlignOptions::getOptions()
{
  if(instance==nullptr) {
    cerr << "AlignOptions::getInstance - Options haven't been initialized"<<endl;
    throw "AlignOptions::getInstance - Options haven't been initialized";
  }

  return instance;
}

AlignOptions* AlignOptions::createOptions(int argc, char* argv[] )
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new AlignOptions(argc, argv);
  return instance;
}

AlignOptions* AlignOptions::createOptions()
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new AlignOptions();
  return instance;
}


