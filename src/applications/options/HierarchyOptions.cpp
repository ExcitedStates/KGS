//
// Created by Dominik Budday on 20.12.16.
//

#include "HierarchyOptions.h"
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

HierarchyOptions::HierarchyOptions(){
  initializeVariables();
}

HierarchyOptions::HierarchyOptions(int argc, char* argv[])
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
    if(arg=="--annotation"){                    annotationFile = argv[++i];                         continue; }
    if(arg=="--hbondMethod"){                   hydrogenbondMethod = argv[++i];                     continue; }
    if(arg=="--hbondFile"){                     hydrogenbondFile = argv[++i];                       continue; }
    if(arg=="--extraCovBonds"){                 Util::split( string(argv[++i]),',', extraCovBonds );   continue; }
    if(arg=="--workingDirectory"){              workingDirectory = argv[++i];                       continue; }
    if(arg=="--samples" || arg=="-s"){          samples = atoi(argv[++i]);                          continue; }
    if(arg=="--collisionFactor" || arg=="-c"){  collisionFactor = atof(argv[++i]);                  continue; }
    if(arg=="--stepSize"){                      stepSize = atof(argv[++i]);                         continue; }
    if(arg=="--saveData"){                      saveData = atoi(argv[++i]);                         continue; }
    if(arg=="--residueNetwork" || arg=="-res"){ residueNetwork = argv[++i];                         continue; }
    if(arg=="--root"){                          Util::split( string(argv[++i]),',', roots );        continue; }
    if(arg=="--collisionCheck"){                collisionCheck = argv[++i];                         continue; }
    if(arg=="--svdCutoff"){                     svdCutoff = atof(argv[++i]);                         continue; }

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

  //Set workingDirectory and moleculeName using the initialStructureFile.
  char* tmp = realpath(initialStructureFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<initialStructureFile<<" is not a valid PDB-file"<<endl; exit(-1); }
  string pdb_file(tmp);
  int nameSplit = pdb_file.find_last_of("/\\");
  if(workingDirectory.empty()) {
    workingDirectory = pdb_file.substr(0, nameSplit + 1);
    log("so")<<"Changing working directory to "<<workingDirectory<<" using initial"<<endl;
  }

  free(tmp);
}


void HierarchyOptions::initializeVariables(){
  initialStructureFile      = "";
  annotationFile            = "";
  hydrogenbondFile          = "";
  hydrogenbondMethod        = "";
  workingDirectory          = "";
  samples                   = 0;
  collisionFactor           = 0.75;
  stepSize                  = 1.0;
  saveData                  =  1;
  residueNetwork            = "all";
  roots                     = {1}; //Choose the first atom
  collisionCheck            = "all";
  svdCutoff                 = 1.0e-12;
}

void HierarchyOptions::print(){
  log("so")<<"Sampling options:"<<std::endl;
  log("so")<<"\t--initial "<<initialStructureFile<<endl;
//  log("so")<<"\t--target "<<targetStructureFile<<endl;
  log("so")<<"\t--annotation "<<annotationFile<<endl;
  if(!hydrogenbondFile.empty()) {
    log("so")<<"\t--hbondMethod "<<hydrogenbondMethod<<endl;
    log("so")<<"\t--hbondFile "<<hydrogenbondFile<<endl;
  }
  log("so")<<"\t--extraCovBonds "; for(unsigned int i=0;i<extraCovBonds.size();i++) log("so")<<extraCovBonds[i]<<","; log("so")<<endl;
  log("so")<<"\t--workingDirectory "<<workingDirectory<<endl;
  log("so")<<"\t--samples "<<samples<<endl;
  log("so")<<"\t--collisionFactor "<<collisionFactor<<endl;
  log("so")<<"\t--saveData "<<saveData<<endl;
  log("so")<<"\t--root "; for(unsigned int i=0;i<roots.size();i++) log("so")<<roots[i]<<" "; log("so")<<endl;
  log("so")<<"\t--collisionCheck "<<collisionCheck<<endl;
  log("so")<<"\t--svdCutoff "<<svdCutoff<<endl;
}

void HierarchyOptions::printUsage(char* pname){
  log("so")<<"Usage: "<<pname<<" --initial <pdb-file> [option list]"<<endl;
  log("so")<<"Performs h-bond hierarchy analysis for atom positions and hydrogen bonds in initial file"<<endl;
  log("so")<<endl;
  log("so")<<"Options:"<<endl;

  log("so")<<"  --initial <pdb-file> \t: Specifies the initial structure."<<endl;
  log("so")<<"  --annotation <file-path> \t: Annotations can specify secondary structures or other things ";
//  log("so")<<"\t--hbondMethod <user|dssr|rnaview|first|kinari|hbplus|vadar|identify> \t: Format of the --hbondFile. If no --hbondFile argument is provided, instructions ";
//  log("so")<<"how to generate a hbondFile are printed."<<endl;
//  log("so")<<"\t--hbondFile <path to hydrogen bond file> \t: Hydrogen bond definition file. The format is specified by the choice ";
//  log("so")<<"of --hbondMethod. Leave this field blank for instructions how to generate the hbond file."<<endl;
  log("so")<<"  --extraCovBonds <resi1>/<name1>-<resi2>/<name2>[,...] \t: Extra covalent bonds. Can override an h-bond."<<endl;
  log("so")<<"  --workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;
  log("so")<<"  --samples <int> \t: Samples to generate: 0 (default, no samples) up to # of dofs "<<endl;
  log("so")<<"  --collisionFactor, -c <real number> \t: A number that is multiplied with the van der Waals radius when ";
  log("so")<<"checking for collisions. The default is 0.75."<<endl;
  log("so")<<"  --saveData <0|1|2|3>\t: Indicate whether files shall be saved! 0=none, 1=pdb, 2=pdb and q, 3=all. Default: 1"<<endl;
  log("so")<<"  --roots <int>[,<int>..]\t: Atom IDs of chain roots. Defaults to first atom of each chain."<<endl;
  log("so")<<"  --collisionCheck <string>\t: Atoms used for collision detection: all (default), heavy, backbone"<<endl;
  log("so")<<"  --svdCutoff <real number> \t: Smallest singular value considered as part of the nullspace. Default: 1.0e-12. Higher value can artificially increase nullspace."<<endl;
}



/** From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
inline bool HierarchyOptions::fileExists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

HierarchyOptions* HierarchyOptions::instance = nullptr;

HierarchyOptions* HierarchyOptions::getOptions()
{
  if(instance==nullptr) {
    cerr << "HierarchyOptions::getInstance - Sampling options haven't been initialized"<<endl;
    throw "HierarchyOptions::getInstance - Sampling options haven't been initialized";
  }

  return instance;
}

HierarchyOptions* HierarchyOptions::createOptions(int argc, char* argv[] )
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new HierarchyOptions(argc, argv);
  return instance;
}

HierarchyOptions* HierarchyOptions::createOptions()
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new HierarchyOptions();
  return instance;
}
