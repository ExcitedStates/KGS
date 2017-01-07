//
// Created by Dominik Budday on 19.12.16.
//

#include "RigidityOptions.h"
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

RigidityOptions::RigidityOptions(){
  initializeVariables();
}

RigidityOptions::RigidityOptions(int argc, char* argv[])
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
    if(arg=="--collisionFactor" || arg=="-c"){  collisionFactor = atof(argv[++i]);                  continue; }
    if(arg=="--seed"){                          seed = atoi(argv[++i]);                             continue; }
    if(arg=="--saveData"){                      saveData = atoi(argv[++i]);                         continue; }
    if(arg=="--residueNetwork" || arg=="-res"){ residueNetwork = argv[++i];                         continue; }
    if(arg=="--preventClashes"){                preventClashes = Util::stob(argv[++i]);             continue; }
    if(arg=="--root"){                          Util::split( string(argv[++i]),',', roots );        continue; }
    if(arg=="--collisionCheck"){                collisionCheck = argv[++i];                         continue; }
    if(arg=="--svdCutoff"){                     svdCutoff = atof(argv[++i]);                         continue; }
    if(arg=="--collapseRigidEdges"){             collapseRigid = atoi(argv[++i]);                    continue; }
//    if(arg=="--relativeDistances"){             relativeDistances = argv[++i];                      continue; }

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

  if(collapseRigid<0 || collapseRigid>2){
    log("so")<<endl<<"--collapseRigidEdges must be an integer between 0 and 2 (is "<<collapseRigid<<")"<<endl;
  }

  free(tmp);
}


void RigidityOptions::initializeVariables(){
  initialStructureFile      = "";
//  targetStructureFile       = "";
  annotationFile            = "";
  hydrogenbondFile          = "";
  hydrogenbondMethod        = "";
  workingDirectory          = "";
  collisionFactor           = 0.75;
  seed                      =  418;
  saveData                  =  1;
  preventClashes            = false; //ToDo: include steric constraints in rigidity analysis
  residueNetwork            = "all";
  roots                     = {1}; //Choose the first atom
  collisionCheck            = "all";
  svdCutoff                 = 1.0e-12;
  collapseRigid             = true;
}

void RigidityOptions::print(){
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
  log("so")<<"\t--collisionFactor "<<collisionFactor<<endl;
  log("so")<<"\t--seed "<<seed<<endl;
  log("so")<<"\t--saveData "<<saveData<<endl;
  log("so")<<"\t--preventClashes "<<preventClashes<<endl;
  log("so")<<"\t--root "; for(unsigned int i=0;i<roots.size();i++) log("so")<<roots[i]<<" "; log("so")<<endl;
  log("so")<<"\t--collisionCheck "<<collisionCheck<<endl;
  log("so")<<"\t--svdCutoff "<<svdCutoff<<endl;
  log("so")<<"\t--collapseRigidEdges "<<collapseRigid<<endl;
}

void RigidityOptions::printUsage(char* pname){
  log("so")<<"Usage: "<<pname<<" --initial <pdb-file> [option list]"<<endl;
  log("so")<<"Performs rigidity analysis for atom positions and hydrogen bonds in initial file"<<endl;
  log("so")<<endl;
  log("so")<<"Options:"<<endl;

  log("so")<<"\t--initial <pdb-file> \t: Specifies the initial structure."<<endl;

  log("so")<<"\t--annotation <file-path> \t: Annotations can specify secondary structures or other things ";
  log("so")<<"relevant to the sampling strategy. For RNA, standard WC will indicate non-free residues that wont be rebuilt"<<endl;

//  log("so")<<"\t--hbondMethod <user|dssr|rnaview|first|kinari|hbplus|vadar|identify> \t: Format of the --hbondFile. If no --hbondFile argument is provided, instructions ";
//  log("so")<<"how to generate a hbondFile are printed."<<endl;

//  log("so")<<"\t--hbondFile <path to hydrogen bond file> \t: Hydrogen bond definition file. The format is specified by the choice ";
//  log("so")<<"of --hbondMethod. Leave this field blank for instructions how to generate the hbond file."<<endl;
  log("so")<<"\t--extraCovBonds <resi1>/<name1>-<resi2>/<name2>[,...] \t: Extra covalent bonds. Can override an h-bond."<<endl;

  log("so")<<"\t--workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;

  log("so")<<"\t--collisionFactor, -c <real number> \t: A number that is multiplied with the van der Waals radius when ";
  log("so")<<"checking for collisions. The default is 0.75."<<endl;

  log("so")<<"\t--seed <integer>\t: The seed used for random number generation (Standard: 418)"<<endl;

  log("so")<<"\t--saveData <0|1|2>\t: Indicate whether files shall be saved! 0=none, 1=pdb and q, 2=all. Default: 1"<<endl;

  log("so")<<"\t--preventClashes <true|false>\t: Use clashing atoms to define additional constraints and prevent clash. Default false (currently not enabled for rigidity)."<<endl;

  log("so")<<"\t--residueNetwork <selection-pattern>\t: A pymol-like pattern that specifies mobile residues in rigidity analysis. Default is 'all', other options currently not enabled."<<endl;

  log("so")<<"\t--roots <comma-sep list of int>\t: The atom ID which will be part of the root rigid bodies. Specify one for each chain, as comma-separated list of ints."<<endl;

  log("so")<<"\t--collisionCheck <string>\t: atoms used for collision detection: all (default), heavy, backbone"<<endl;
//
  log("so")<<"\t--svdCutoff <real number> \t: Smallest singular value considered as part of the nullspace, default 1.0e-12. Higher value can artificially increase nullspace."<<endl;

  log("so")<<"\t--collapseRigidEdges <0|1|2> \t: Indicates whether to speed up null-space computation by collapsing rigid edges. 0: Dont collapse. 1: Collapse covalent bonds. 2: Collapse covalent and hydrogen bonds. Default 2. (other options currently not enabled)."<<endl;
}



/** From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
inline bool RigidityOptions::fileExists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

RigidityOptions* RigidityOptions::instance = nullptr;

RigidityOptions* RigidityOptions::getOptions()
{
  if(instance==nullptr) {
    cerr << "RigidityOptions::getInstance - Sampling options haven't been initialized"<<endl;
    throw "RigidityOptions::getInstance - Sampling options haven't been initialized";
  }

  return instance;
}

RigidityOptions* RigidityOptions::createOptions(int argc, char* argv[] )
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new RigidityOptions(argc, argv);
  return instance;
}

RigidityOptions* RigidityOptions::createOptions()
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new RigidityOptions();
  return instance;
}