
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
#include <moves/ClashAvoidingMove.h>
#include <metrics/Dihedral.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <moves/DecreaseStepMove.h>
#include <loopclosure/ExactIK.h>
#include <metrics/RMSDnosuper.h>

using namespace std;

extern double jacobianTime;
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

  //SamplingOptions options(argc,argv);
  SamplingOptions::createOptions(argc, argv);

  SamplingOptions &options = *(SamplingOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//SamplingOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  string pdb_file = options.initialStructureFile;
  Molecule protein;
  protein.setCollisionFactor(options.collisionFactor);

  IO::readPdb( &protein, pdb_file, options.extraCovBonds );
  options.setResidueNetwork(&protein);
  IO::readHbonds( &protein, options.hydrogenbondFile );
  IO::readRigidbody( &protein, options.residueNetwork );
  protein.buildSpanningTree();//with the rigid body tree in place, we can generate a configuration
  protein.setConfiguration(new Configuration(&protein));
  protein.m_initialCollisions = protein.getAllCollisions();

  if(options.selectionMoving.empty()){
    cerr<<"Must supply --selectionMoving (e.g. \"resi 17-22 and resi 50-55\")"<<endl;
    exit(-1);
  }


  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein.getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein.m_initialCollisions.size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->CycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein.m_spanning_tree->getNumDOFs() << " DOFs of which " << protein.m_spanning_tree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  //Initialize metric
  metrics::Metric* metric = nullptr;
  try {
    Selection metricSelection(options.metricPattern);
    if(SamplingOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move = nullptr;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new ClashAvoidingMove();
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(SamplingOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  //Initialize direction
  Direction* direction = nullptr;
  if(options.gradient == 0)      direction = new RandomDirection();
  else if(options.gradient <= 2) direction = new DihedralDirection();
  else if(options.gradient <= 4) direction = new MSDDirection();
  else if(options.gradient <= 5) direction = new LSNullspaceDirection();
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
          Residue* res1 = protein.getChain("A")->getResidue(r1);
          Residue* res2 = protein.getChain("A")->getResidue(r2);
          Residue* res3 = protein.getChain("A")->getResidue(r3);
          if(ik.validRebuildLoop(res1,res2,res3))
            ikTriples.push_back(make_tuple(res1, res2, res3));
        }
      }
    }
  }

  //Initialize planner
  SamplingPlanner* planner = nullptr;
  if(options.planner_string=="binnedrrt")         planner = new RRTPlanner(     &protein, *move, *metric, *direction );
  else if(options.planner_string.empty())         planner = new RRTPlanner(     &protein, *move, *metric, *direction );
  else if(options.planner_string=="dihedralrrt")  planner = new DihedralRRT(    &protein, *move, *metric, *direction );
  else if(options.planner_string=="poisson")      planner = new PoissonPlanner( &protein, *move, *metric );
  else if(options.planner_string=="poisson2")     planner = new PoissonPlanner2(&protein, *move, *metric, ikTriples );
  else{
    cerr<<"Unknown --planner option specified!"<<endl;
    exit(-1);
  }

  if(options.saveData > 1){
    string out = options.workingDirectory + "output/" + protein.getName() + "_q_0.txt";
    IO::writeQ(&protein, protein.m_conf, out);
  }

  log("samplingStatus")<<"Sampling ...\n";
  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();

  //Start exploring
  planner->GenerateSamples();

  //Print final status
  double end_time = timer.ElapsedTime();
  std::list<Configuration*> m_samples = planner->Samples();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(m_samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianTime<<" seconds\n";
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

