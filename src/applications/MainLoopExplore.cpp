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



#include <stdexcept>
#include <string>
#include <iostream>
#include <list>
#include "core/Configuration.h"
#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"

#include "planners/SamplingPlanner.h"
#include "planners/RRTPlanner.h"
#include "planners/DihedralRRT.h"
#include "planners/PoissonPlanner.h"
#include "planners/PoissonPlanner2.h"
#include "planners/BidirectionalMovingFront.h"
#include <moves/NullspaceMove.h>
#include <moves/FastClashAvoidingMove.h>
#include <metrics/Dihedral.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <moves/DecreaseStepMove.h>
#include <loopclosure/ExactIK.h>
#include <metrics/RMSDnosuper.h>
#include <applications/options/ExploreOptions.h>

using namespace std;

extern double jacobianAndNullspaceTime;
extern double rigidityTime;
extern double selectNodeTime;

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

  ofstream reportStream;
  reportStream.open("kgs_report.log");
  enableLogger("report", reportStream);

  ofstream debugStream;
  debugStream.open("kgs_debug.log");
  enableLogger("debug", debugStream);

  ofstream plannerStream;
  plannerStream.open("kgs_planner.log");
  enableLogger("dominik", plannerStream);

  //ExploreOptions options(argc,argv);
  ExploreOptions::createOptions(argc, argv);

  ExploreOptions &options = *(ExploreOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//ExploreOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  Selection resNetwork(options.gradientSelection);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      options.extraCovBonds,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  protein->initializeTree(resNetwork,options.collisionFactor,options.roots);

  if(options.gradientSelection.empty()){
    cerr<<"Must supply --gradientSelection (e.g. \"resi 17-22 and resi 50-55\")"<<endl;
    exit(-1);
  }

//  if(options.hydrogenbondMethod=="user")
//    IO::readHbonds( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="rnaview")
//    IO::readHbonds_rnaview( protein, options.hydrogenbondFile, options.annotationFile.empty() );
//  else if(options.hydrogenbondMethod=="first" || options.hydrogenbondMethod=="FIRST")
//    IO::readHbonds_first( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="kinari" || options.hydrogenbondMethod=="KINARI")
//    IO::readHbonds_kinari( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="hbplus" || options.hydrogenbondMethod=="hbPlus")
//    IO::readHbonds_hbPlus( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="vadar")
//    IO::readHbonds_vadar( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="dssr")
//    IO::readHbonds_dssr( protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="identify")
//    HbondIdentifier::identifyHbonds(protein);


  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein->getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  try {
    Selection metricSelection(options.metricSelection);
    if(ExploreOptions::getOptions()->metric_string=="rmsd")        metric = new metrics::RMSD(metricSelection);
    if(ExploreOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(ExploreOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move = nullptr;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new FastClashAvoidingMove( options.maxRotation,
                                      options.decreaseSteps,
                                      options.collisionCheck,
                                      options.projectConstraints );
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(ExploreOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  //Initialize m_direction
  Direction* direction = nullptr;
  if(options.gradient == 0)      direction = new RandomDirection(resNetwork);
  else if(options.gradient <= 2) direction = new DihedralDirection(resNetwork);
  else if(options.gradient <= 4) direction = new MSDDirection(resNetwork);
  else if(options.gradient <= 5) direction = new LSNullspaceDirection(resNetwork);
  else{
    cerr<<"Invalid --gradient specified. Must be between 0 and 5"<<endl;
    exit(-1);
  }

  vector<ResTriple> ikTriples;
  //vector<pair<int, int> > intervals = {{20,  24}};
  vector<pair<int, int> > intervals = {{14,  24}};
//  vector<pair<int, int> > intervals = {{20,  24},
//                                       {116, 121},
//                                       {124, 127}};

  //Generate all triples within intervals
  ExactIK ik;
  for( auto ival: intervals ) {
    for( int r1=ival.first; r1<=ival.second; r1++ ) {
      for( int r2=r1 + 1; r2<=ival.second; r2++ ) {
        for( int r3=r2 + 1; r3<=ival.second; r3++ ) {
          Residue* res1 = protein->getChain("A")->getResidue(r1);
          Residue* res2 = protein->getChain("A")->getResidue(r2);
          Residue* res3 = protein->getChain("A")->getResidue(r3);
          if(ik.validRebuildLoop(res1,res2,res3))
            ikTriples.push_back(make_tuple(res1, res2, res3));
        }
      }
    }
  }

  //Initialize planner
  SamplingPlanner* planner = nullptr;
  if(options.planner_string.empty() || options.planner_string=="binnedrrt")
    planner = new RRTPlanner(
        protein,
        direction,
        options.explorationRadius,
        options.samplesToGenerate,
        options.gradient,
        options.stepSize,
        options.maxRotation,
        options.scaleToRadius
    );
  else if(options.planner_string=="dihedralrrt")
    planner = new DihedralRRT(
        protein,
        direction,
        options.samplesToGenerate,
        options.explorationRadius,
        options.maxRotation,
        options.sampleRandom
    );
  else if(options.planner_string=="poisson")
    planner = new PoissonPlanner(
        protein,
        options.samplesToGenerate,
        options.poissonMaxRejectsBeforeClose,
        options.stepSize,
        options.gradientSelection
    );
  else if(options.planner_string=="poisson2")
    planner = new PoissonPlanner2(
        protein,
        ikTriples,
        options.samplesToGenerate,
        options.poissonMaxRejectsBeforeClose,
        options.stepSize,
        options.residueNetwork
    );
  else{
    cerr<<"Unknown --planner option specified!"<<endl;
    exit(-1);
  }
  planner->initialize(move, metric, options.workingDirectory, options.saveData);

  if(options.saveData > 1){
    string out = options.workingDirectory + "output/" + protein->getName() + "_q_0.txt";
    IO::writeQ(protein, protein->m_conf, out);
  }

  log("samplingStatus")<<"Sampling ...\n";
  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  //Start exploring
  planner->generateSamples();

  //Print final status
  double end_time = timer.ElapsedTime();
  std::list<Configuration*> m_samples = planner->getSamples();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianAndNullspaceTime<<" seconds\n";
  log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
  log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";


  if(options.saveData > 0){
    planner->createTrajectory();
  }

  reportStream.close();
  plannerStream.close();

  //Clean up
  delete planner;
  delete direction;

  return 0;
}

