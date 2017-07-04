#include "RelativeTransitionOptions.h"
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

RelativeTransitionOptions::RelativeTransitionOptions():
    ApplicationOptions()
{
  initializeVariables();
}

RelativeTransitionOptions::RelativeTransitionOptions(int argc, char* argv[]):
    ApplicationOptions(argc, argv)
{

  initializeVariables();

  if(argc<2){
    enableLogger("so");
    printUsage(argv[0]);
    exit(-1);
  }

  for(int i=1; i < argc; i++){
    string arg = argv[i];
    if(arg=="--initial"){                       initialStructureFile = argv[++i];                    continue; }
    if(arg=="--workingDirectory"){              workingDirectory = argv[++i];                        continue; }
    if(arg=="--samples" || arg=="-s"){          samplesToGenerate = atoi(argv[++i]);                 continue; }
    if(arg=="--hbondMethod"){                   hydrogenbondMethod = argv[++i];                      continue; }
    if(arg=="--hbondFile"){                     hydrogenbondFile = argv[++i];                        continue; }
    if(arg=="--extraCovBonds"){                 Util::split( string(argv[++i]),',', extraCovBonds ); continue; }
    if(arg=="--gradient" || arg=="-g"){         gradient = atoi(argv[++i]);                          continue; }
    if(arg=="--collisionFactor" || arg=="-c"){  collisionFactor = atof(argv[++i]);                   continue; }
    if(arg=="--decreaseSteps" ){                decreaseSteps = atoi(argv[++i]);                     continue; }
    if(arg=="--decreaseFactor"){                decreaseFactor = atof(argv[++i]);                    continue; }
    if(arg=="--stepSize"){                      stepSize = atof(argv[++i]);                          continue; }
    if(arg=="--maxRotation"   ){                maxRotation = atof(argv[++i]);                       continue; }
    if(arg=="--metric"){                        metric_string = argv[++i];                           continue; }
    if(arg=="--metricSelection"){               metricSelection = argv[++i];                         continue; }
    if(arg=="--seed"){                          seed = atoi(argv[++i]);                              continue; }
    if(arg=="--biasToTarget" || arg=="-bias"){  biasToTarget = atof(argv[++i]);                      continue; }
    if(arg=="--convergeDistance"){              convergeDistance = atof(argv[++i]);                  continue; }
    if(arg=="--residueNetwork" || arg=="-res"){ residueNetwork = argv[++i];                          continue; }
    if(arg=="--preventClashes"){                preventClashes = Util::stob(argv[++i]);              continue; }
    if(arg=="--alignSelection"){                alignSelection = argv[++i];                          continue; }
    if(arg=="--gradientSelection"){             gradientSelection = argv[++i];                       continue; }
    if(arg=="--roots"){                         Util::split( string(argv[++i]),',', roots );         continue; }
    if(arg=="--projectConstraints"){            projectConstraints = Util::stob(argv[++i]);          continue; }
    if(arg=="--collisionCheck"){                collisionCheck = argv[++i];                          continue; }
    if(arg=="--svdCutoff"){                     svdCutoff = atof(argv[++i]);                         continue; }
    if(arg=="--collapseRigidEdges"){            collapseRigid = atoi(argv[++i]);                     continue; }
    if(arg=="--relativeDistances"){             relativeDistances = argv[++i];                       continue; }
    if(arg=="--predictStrain"){                 predictStrain = Util::stob(argv[++i]);               continue; }
    if(arg=="--logger"){                        ++i;                                                 continue; }

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


  //Set workingDirectory and moleculeName using the initialStructureFile.
  char* tmp = realpath(initialStructureFile.c_str(), nullptr);
  if(tmp==nullptr){ cerr<<initialStructureFile<<" is not a valid PDB-file"<<endl; exit(-1); }
  string pdb_file(tmp);
  int nameSplit = pdb_file.find_last_of("/\\");
  if(workingDirectory.empty()) {
    workingDirectory = pdb_file.substr(0, nameSplit + 1);
    log("so")<<"Changing working directory to "<<workingDirectory<<" using initial"<<endl;
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

  //Ensure that decreaseSteps>0 if preventClashes set
  if(preventClashes && decreaseSteps==0)
    decreaseSteps=5;

  if(collapseRigid<0 || collapseRigid>2){
    log("so")<<endl<<"--collapseRigidEdges must be an integer between 0 and 2 (is "<<collapseRigid<<")"<<endl;
  }


  free(tmp);
}


void RelativeTransitionOptions::initializeVariables(){
  initialStructureFile      = "";
  workingDirectory          = "";
  samplesToGenerate         = 10;
  hydrogenbondFile          = "";
  hydrogenbondMethod        = "";
  gradient                  = 5;
  collisionFactor           = 0.75;
  decreaseSteps             = 0;
  decreaseFactor            = 0.5;
  stepSize                  = 1;
  maxRotation               = 3.1415/18;
  metric_string             = "rmsd";
  metricSelection           = "heavy";
  seed                      =  418;
  biasToTarget              = 0.5;
  convergeDistance          = -1.0; ///<Changes depending on m_metric. Initialize <0, it is set depending on metric
  preventClashes            = false;
  alignSelection            = "heavy";
  gradientSelection         = "heavy";
  residueNetwork            = "all";
  roots                     = {1}; //Choose the first atom
  projectConstraints        = true;
  collisionCheck            = "all";
  svdCutoff                 = 1.0e-12;
  collapseRigid             = false;
  relativeDistances         = "";
  predictStrain             = false;
}

void RelativeTransitionOptions::print(){
  log("so")<<"Sampling options:"<<std::endl;
  log("so")<<"  --initial "<<initialStructureFile<<endl;
  log("so")<<"  --workingDirectory "<<workingDirectory<<endl;
  if(!hydrogenbondFile.empty()) {
    log("so")<<"  --hbondMethod "<<hydrogenbondMethod<<endl;
    log("so")<<"  --hbondFile "<<hydrogenbondFile<<endl;
  }
  log("so")<<"  --samples "<<samplesToGenerate<<endl;
  log("so")<<"  --extraCovBonds "; for(unsigned int i=0;i<extraCovBonds.size();i++) log("so")<<extraCovBonds[i]<<","; log("so")<<endl;
  log("so")<<"  --gradient "<<gradient<<endl;
  log("so")<<"  --collisionFactor "<<collisionFactor<<endl;
  log("so")<<"  --decreaseSteps "<< decreaseSteps <<endl;
  log("so")<<"  --decreaseFactor "<<decreaseFactor<<endl;
  log("so")<<"  --stepSize "<<stepSize<<endl;
  log("so")<<"  --maxRotation "<<maxRotation<<endl;
  log("so")<<"  --metric "<<metric_string<<endl;
  log("so")<<"  --metricSelection "<<metricSelection<<endl;
  log("so")<<"  --seed "<<seed<<endl;
  log("so")<<"  --biasToTarget "<<biasToTarget<<endl;
  log("so")<<"  --convergeDistance "<<convergeDistance<<endl;
  log("so")<<"  --preventClashes "<<preventClashes<<endl;
  log("so")<<"  --gradientSelection "<<gradientSelection<<endl;
  log("so")<<"  --roots "; for(unsigned int i=0;i<roots.size();i++) log("so")<<roots[i]<<" "; log("so")<<endl;
  log("so")<<"  --collisionCheck "<<collisionCheck<<endl;
  log("so")<<"  --svdCutoff "<<svdCutoff<<endl;
  log("so")<<"  --collapseRigidEdges "<<collapseRigid<<endl;
  log("so")<<"  --relativeDistances "<<relativeDistances<<endl;
  log("so")<<"  --predictStrain "<<predictStrain<<endl;
}

void RelativeTransitionOptions::printUsage(char* pname){
  log("so")<<"Usage: "<<pname<<" explore [option list]"<<endl;
  log("so")<<"The KGS program will start sampling using the specified options."<<endl;
  log("so")<<"Options are:"<<endl;
  log("so")<<"  --initial <path to structure> \t: Specifies the initial structure."<<endl;
  log("so")<<"  --workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;
  log("so")<<"  --extraCovBonds <int>-<int>[,...] \t: Extra covalent bonds specified with atom-IDs. Can override an h-bond."<<endl;
  log("so")<<"  --samples, -s <whole number> \t: Indicates the number of samples to generate. Default is 10."<<endl;
  log("so")<<"  --gradient <0|1|2|3|4\5> \t: Indicates method to calculate a new perturbation or gradient. 0 = random; 1 = torsion, no blending; 2 = torsion, low pass blending; 3 = msd, no blending; 4 = msd, low pass blending; 5 = least square gradient. Default is 0."<<endl;
  log("so")<<"  --collisionFactor, -c <real number> \t: A number that is multiplied with the van der Waals radius when checking for collisions. The default is 0.75."<<endl;
  log("so")<<"  --stepSize <real number> \t: Initial step size to next sample as the norm of dihedral perturbation. Default 1."<<endl;
  log("so")<<"  --maxRotation <real number> \t: The largest allowable change in torsion angle. The default is 10°."<<endl;
  log("so")<<"  --metric <rmsd|rmsdnosuper|dihedral> \t: The metric to use in sampler. Default is 'rmsd'."<<endl;
  log("so")<<"  --metricSelection <selection-pattern>\t: A pymol-like pattern that indicates which subset of atoms the metric operates on. Default is 'heavy'."<<endl;
  log("so")<<"  --seed <integer>\t: The seed used for random number generation (Standard: 418)"<<endl;
  log("so")<<"  --biasToTarget, -bias <real number> \t: Percentage mixing the directed move with a random move. Default 0.5."<<endl;
  log("so")<<"  --convergeDistance <real number> \t: The distance under which the goal conformation is considered reached. Default is 0.1 for RMSD and 1e-8 for Dihedral metric."<<endl;
  log("so")<<"  --preventClashes\t: Use clashing atoms to define additional constraints and prevent clash. Default: true."<<endl;
  log("so")<<"  --gradientSelection <selection-pattern>\t: A pymol-like pattern that pecifies the residues of the molecule that are used to determine the gradient. Default is <heavy>."<<endl;
  log("so")<<"  --residueNetwork <selection-pattern>\t: A pymol-like pattern that specifies mobile residues during sampling (e.g. limited to single flexible loop). Default is 'all'."<<endl;
  log("so")<<"  --roots <comma-sep list of int>\t: The atom ID which will be part of the root rigid bodies. Specify one for each chain, as comma-separated list of ints."<<endl;
  log("so")<<"  --collisionCheck <string>\t: atoms used for collision detection: all (default), heavy, backbone"<<endl;
  log("so")<<"  --svdCutoff <real number> \t: Smallest singular value considered as part of the nullspace, default 1.0e-12."<<endl;
  log("so")<<"  --collapseRigidEdges <0|1|2> \t: Indicates whether to speed up null-space computation by collapsing rigid edges. 0: Dont collapse. 1: Collapse covalent bonds. 2: Collapse covalent and hydrogen bonds. Default 0."<<endl;
  log("so")<<"  --relativeDistances <file name> \t: File with the desired distance between atoms of residueNetwork option. Each line is one couple 'id atom_id2,id atom_id2,distance'. Each line should contain only two spaces, one after each word 'id'"<<endl;
  log("so")<<"  --predictStrain <true|false> \t: Indicate whether strain on constraints should be printed."<<endl;
}

