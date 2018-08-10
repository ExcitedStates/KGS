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

//
// Created by Dominik Budday on 02.08.18.
//

#include "DeerOptions.h"

#include "DeerOptions.h"
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
#include "ApplicationOptions.h"

using namespace std;

DeerOptions::DeerOptions():
    ApplicationOptions()
{
  initializeVariables();
}

DeerOptions::DeerOptions(int argc, char* argv[]):
    ApplicationOptions(argc, argv)
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
    if(arg=="--hbondMethod"){                   hydrogenbondMethod = argv[++i];                     continue; }
    if(arg=="--hbondFile"){                     hydrogenbondFile = argv[++i];                       continue; }
    if(arg=="--extraCovBonds"){                 Util::split( string(argv[++i]),',', extraCovBonds );   continue; }
    if(arg=="--workingDirectory"){              workingDirectory = argv[++i];                       continue; }
    if(arg=="--samples" || arg=="-s"){          samplesToGenerate = atoi(argv[++i]);                continue; }
    if(arg=="--gradient" || arg=="-g"){         gradient = atoi(argv[++i]);                         continue; }
    if(arg=="--collisionFactor" || arg=="-c"){  collisionFactor = atof(argv[++i]);                  continue; }
    if(arg=="--decreaseSteps" ){                decreaseSteps = atoi(argv[++i]);                    continue; }
    if(arg=="--decreaseFactor"){                decreaseFactor = atof(argv[++i]);                   continue; }
    if(arg=="--stepSize"){                      stepSize = atof(argv[++i]);                         continue; }
    if(arg=="--maxRotation"   ){                maxRotation = atof(argv[++i]);                      continue; }
    if(arg=="--metric"){                        metric_string = argv[++i];                          continue; }
    if(arg=="--metricSelection"){               metricSelection = argv[++i];                        continue; }
    if(arg=="--planner"   ){                    planner_string = argv[++i];                         continue; }
    if(arg=="--seed"){                          seed = atoi(argv[++i]);                             continue; }
    if(arg=="--biasToTarget" || arg=="-bias"){  biasToTarget = atof(argv[++i]);                     continue; }
    if(arg=="--convergeDistance"){              convergeDistance = atof(argv[++i]);                 continue; } //Previously rmsdThreshold
    if(arg=="--saveData"){                      saveData = atoi(argv[++i]);                         continue; }
    if(arg=="--residueNetwork" || arg=="-res"){ residueNetwork = argv[++i];                         continue; }
    if(arg=="--preventClashes"){                preventClashes = Util::stob(argv[++i]);             continue; }
    if(arg=="--alignSelection"){                alignSelection = argv[++i];                         continue; }
    if(arg=="--gradientSelection"){             gradientSelection = argv[++i];                      continue; }
    if(arg=="--root"){                          Util::split( string(argv[++i]),',', roots );        continue; }
    if(arg=="--projectConstraints"){            projectConstraints = Util::stob(argv[++i]);         continue; }
    if(arg=="--collisionCheck"){                collisionCheck = argv[++i];                         continue; }
    if(arg=="--frontSize"){                     frontSize = atoi(argv[++i]);                        continue; }
    if(arg=="--svdCutoff"){                     svdCutoff = atof(argv[++i]);                        continue; }
    if(arg=="--collapseRigidEdges"){            collapseRigid = atoi(argv[++i]);                    continue; }
    if(arg=="--relativeDistances"){             relativeDistances = argv[++i];                      continue; }
    if(arg=="--predictStrain"){                 predictStrain = Util::stob(argv[++i]);              continue; }
    if(arg=="--logger"){                        ++i;                                                continue; }
    if(arg=="--radius" || arg=="-r"){           explorationRadius = atof(argv[++i]);                continue; }
    if(arg=="--explore"){                       explore = Util::stob(argv[++i]);                    continue; }

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
    cerr<<"No initial structure file supplied"<<endl<<endl;
    exit(-1);
  }

  //Check relativeDistances
  if(relativeDistances.empty()) {
    enableLogger("so");
    cerr<<endl<<"No relative distance file supplied"<<endl<<endl;
    exit(-1);
  }

  //Check for valid metric
  std::transform(metric_string.begin(), metric_string.end(), metric_string.begin(), ::tolower);
  if(metric_string!="dihedral" && metric_string!="rmsd" && metric_string!="rmsdnosuper"){
    cerr<<"Invalid --metric option: "<<metric_string<<". Must be either 'dihedral', 'RMSD', or 'RMSDnosuper."<<endl;
    exit(-1);
  }

  //Set default convergeDistance
  if(convergeDistance<0.0) {
    if(metric_string=="rmsd")        convergeDistance = 0.01;
    if(metric_string=="rmsdnosuper") convergeDistance = 0.01;
    if(metric_string=="dihedral")    convergeDistance = 1.0e-8;
  }

  //Check for valid planner
  std::transform(planner_string.begin(), planner_string.end(), planner_string.begin(), ::tolower);
  if( planner_string!="dihedralrrt" &&
      planner_string!="binnedrrt" &&
      planner_string!="dccrrt" &&
      planner_string!="deer" &&
      planner_string!="poisson2" &&
      planner_string!="poisson"){
    cerr<<"Invalid --planner option: "<<planner_string<<endl;
    exit(-1);
  }

  //Ensure that decreaseSteps>0 if preventClashes set
  if(preventClashes && decreaseSteps==0)
    decreaseSteps=5;

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


void DeerOptions::initializeVariables(){
  initialStructureFile      = "";
  hydrogenbondFile          = "";
  hydrogenbondMethod        = "";
  workingDirectory          = "";
  samplesToGenerate         = 10;
  gradient                  = 1;
  collisionFactor           = 0.75;
  decreaseSteps             = 2;
  decreaseFactor            = 0.5;
  stepSize                  = 1.0;
  maxRotation               = 3.1415/18;
  metric_string             = "dihedral";
  metricSelection           = "heavy";
  planner_string            = "deer";
  seed                      =  418;
  saveData                  =  1;
  biasToTarget              = 0.4;
  convergeDistance          = 1.0; ///<Changes depending on m_metric. Initialize <0, it is set depending on metric
  preventClashes            = true;
  alignSelection            = "heavy";
  gradientSelection         = "heavy";
  residueNetwork            = "all";
  roots                     = {1}; //Choose the first atom
  projectConstraints        = true;
  collisionCheck            = "all";
  frontSize                 = 50;
  svdCutoff                 = 1.0e-12;
  collapseRigid             = false;
  relativeDistances         = "";
  predictStrain             = false;
  explorationRadius         = 2.0;
  explore                   = false;
}

void DeerOptions::print(){
  log("so")<<"Sampling options:"<<std::endl;
  log("so")<<"\t--initial "<<initialStructureFile<<endl;
  if(!hydrogenbondFile.empty()) {
    log("so")<<"\t--hbondMethod "<<hydrogenbondMethod<<endl;
    log("so")<<"\t--hbondFile "<<hydrogenbondFile<<endl;
  }
  log("so")<<"\t--extraCovBonds "; for(unsigned int i=0;i<extraCovBonds.size();i++) log("so")<<extraCovBonds[i]<<","; log("so")<<endl;
  log("so")<<"\t--workingDirectory "<<workingDirectory<<endl;
  log("so")<<"\t--samples "<<samplesToGenerate<<endl;
  log("so")<<"\t--gradient "<<gradient<<endl;
  log("so")<<"\t--collisionFactor "<<collisionFactor<<endl;
  log("so")<<"\t--decreaseSteps "<< decreaseSteps <<endl;
  log("so")<<"\t--decreaseFactor "<<decreaseFactor<<endl;
  log("so")<<"\t--stepSize "<<stepSize<<endl;
  log("so")<<"\t--maxRotation "<<maxRotation<<endl;
  log("so")<<"\t--metric "<<metric_string<<endl;
  log("so")<<"\t--metricSelection "<<metricSelection<<endl;
  log("so")<<"\t--planner "<<planner_string<<endl;
  log("so")<<"\t--seed "<<seed<<endl;
  log("so")<<"\t--biasToTarget "<<biasToTarget<<endl;
  log("so")<<"\t--convergeDistance "<<convergeDistance<<endl;
  log("so")<<"\t--saveData "<<saveData<<endl;
  log("so")<<"\t--preventClashes "<<preventClashes<<endl;
  log("so")<<"\t--alignSelection "<<alignSelection<<endl;
  log("so")<<"\t--gradientSelection "<<gradientSelection<<endl;
  log("so")<<"\t--root "; for(unsigned int i=0;i<roots.size();i++) log("so")<<roots[i]<<" "; log("so")<<endl;
  log("so")<<"\t--projectConstraints "<<projectConstraints<<endl;
  log("so")<<"\t--collisionCheck "<<collisionCheck<<endl;
  log("so")<<"\t--frontSize "<<frontSize<<endl;
  log("so")<<"\t--svdCutoff "<<svdCutoff<<endl;
  log("so")<<"\t--collapseRigidEdges "<<collapseRigid<<endl;
  log("so")<<"\t--relativeDistances "<<relativeDistances<<endl;
  log("so")<<"\t--predictStrain "<<predictStrain<<endl;
  log("so")<<"\t--radius "<<explorationRadius<<endl;
  log("so")<<"\t--explore "<<explore<<endl;
}

void DeerOptions::printUsage(char* pname){
  log("so")<<"Usage: "<<pname<<" --initial <pdb-file> [option list]"<<endl;
  log("so")<<"Creates a transition towards the specified deer sample distance."<<endl;
  log("so")<<endl;
  log("so")<<"Options:"<<endl;

  log("so")<<"  --initial <pdb-file> \t: Specifies the initial structure."<<endl;
  log("so")<<"  --extraCovBonds <atom ID>-<atom ID>[,...] \t: List of extra covalent bonds. Can override h-bonds."<<endl;
  log("so")<<"  --workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;
  log("so")<<"  --samples, -s <whole number> \t: Indicates the number of samples to generate. Default is 10."<<endl;
  log("so")<<"  --gradient <0|1|2|3|4|5|6> \t: Indicates method to calculate a new perturbation or gradient. 0 = random; 1 = torsion, no blending; 2 = torsion, low pass blending; 3 = msd, no blending; 4 = msd, low pass blending. Default is 0."<<endl;
  log("so")<<"  --collisionFactor, -c <real number> \t: A number that is multiplied with the van der Waals radius when ";
  log("so")<<"checking for collisions. The default is 0.75."<<endl;
  //log("so")<<"  --decreaseSteps <whole number> \t: If a non-colliding structure can not be found, try this many times to ";
  //log("so")<<"decrease the stepsize. Default is 1."<<endl;
  //log("so")<<"  --decreaseFactor <real number> \t: If a non-colliding structure can not be found, decrease the stepsize ";
  //log("so")<<"by this factor. Default is 0.5."<<endl;
  log("so")<<"  --stepSize <real number> \t: Initial step size to next sample as the norm of dihedral perturbation. Default 1."<<endl;
  log("so")<<"  --maxRotation <real number> \t: The largest allowable change in torsion angle. The default is 10°."<<endl;
  log("so")<<"  --metric <rmsd|rmsdnosuper|dihedral> \t: The metric to use in sampler. Default is 'rmsd'."<<endl;
  log("so")<<"  --metricSelection <selection-pattern>\t: A pymol-like pattern that indicates which subset of atoms the metric operates on. Default is 'heavy'."<<endl;
  log("so")<<"  --planner <binnedRRT|dihedralRRT|dccrrt|poisson|deer> \t: The planning strategy used to create samples. Default is dccrrt."<<endl;
  log("so")<<"  --seed <integer>\t: The seed used for random number generation (Standard: 418)"<<endl;
  log("so")<<"  --biasToTarget, -bias <real number> \t: Percentage of using target as 'random seed configuration'. Default 0.1."<<endl;
  log("so")<<"  --convergeDistance <real number> \t: The distance under which the goal conformation is considered reached. Default is 0.1 for RMSD and 1e-8 for Dihedral metric."<<endl;
  log("so")<<"  --saveData <0|1|2>\t: Indicate whether files shall be saved! 0=none, 1=pdb and q, 2=all. Default: 1"<<endl;
  log("so")<<"  --preventClashes "<<(preventClashes?"true":"false: Use clashing atoms to define additional constraints and prevent clash. Default true.")<<endl;
  log("so")<<"  --alignSelection < A pymol-like pattern that indicates which subset of atoms are used during alignment (if specified). Default is 'heavy'."<<endl;
  log("so")<<"  --gradientSelection <selection-pattern>\t: A pymol-like pattern that pecifies the residues of the molecule that are used to determine the gradient. Default is 'heavy'."<<endl;
  log("so")<<"  --residueNetwork <selection-pattern>\t: A pymol-like pattern that specifies mobile residues during sampling (e.g. limited to single flexible loop). Default is 'all'."<<endl;
  log("so")<<"  --roots <comma-sep list of int>\t: Atom IDs of chain roots. Defaults to first atom of each chain."<<endl;
  log("so")<<"  --projectConstraints <true/false>\t: If false, then we don't project moves onto the constraint manifold. Only recommended for testing."<<endl;
  log("so")<<"  --collisionCheck <string>\t: atoms used for collision detection: all (default), heavy, backbone"<<endl;
  log("so")<<"  --frontSize <integer>\t: Size of the propagating front of samples in directed sampling."<<endl;
  log("so")<<"  --svdCutoff <real number> \t: Smallest singular value considered as part of the nullspace, default 1.0e-12."<<endl;
  log("so")<<"  --collapseRigidEdges <0|1|2> \t: Indicates whether to speed up null-space computation by collapsing rigid edges. 0: Dont collapse. 1: Collapse covalent bonds. 2: Collapse covalent and hydrogen bonds. Default 0."<<endl;
  log("so")<<"  --relativeDistances <list of double> \t: has to begin by 'double ' followed by doubles seprated by '+' .It corresponds of the desired distance between atoms of residueNetwork option. "<<endl;
  log("so")<<"  --relativeDistances <file name> \t: File with the desired distance between atoms of residueNetwork option. Each line is one couple 'id atom_id2,id atom_id2,distance'. Each line should contain only two spaces, one after each word 'id'"<<endl;
  log("so")<<"  --predictStrain <true|false> \t: Indicate whether strain on constraints should be printed."<<endl;
  log("so")<<"  --radius <real number> \t: Radius of exploration after DEER distances have been reached."<<endl;
  log("so")<<"  --explore "<<(explore?"true":"false: Explore after DEER distances have been reached. Default false.")<<endl;
}