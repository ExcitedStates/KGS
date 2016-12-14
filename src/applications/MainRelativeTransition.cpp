#include <stdexcept>
#include <string>
#include <iostream>
#include <list>
#include <regex>

#include "core/Configuration.h"
#include "core/Molecule.h"
#include "core/Grid.h"
#include "CTKTimer.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "Logger.h"
#include "moves/NullspaceMove.h"
#include "moves/ClashAvoidingMove.h"
#include "directions/RandomDirection.h"
#include "directions/LSNrelativeDirection.h"
#include "moves/DecreaseStepMove.h"
#include "applications/options/RelativeTransitionOptions.h"

using namespace std;

extern double jacobianTime;
extern double rigidityTime;
extern double selectNodeTime;

void scale_gradient(gsl_vector* gradient, Molecule* mol, double maxRotation);
double dist_to_objective(std::vector< std::tuple<Atom*, Atom*, double> > goal_distances);

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");
//  enableLogger("dominik");


  //RelativeTransitionOptions options(argc,argv);
  RelativeTransitionOptions::createOptions(argc, argv);

  RelativeTransitionOptions &options = *(RelativeTransitionOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//RelativeTransitionOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  string pdb_file = options.initialStructureFile;
  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      pdb_file,
      movingResidues,
      options.extraCovBonds,
      options.roots,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  log() << "Molecule has " << protein->getAtoms().size() << " atoms\n";

  protein->setCollisionFactor(options.collisionFactor);

  if(options.collapseRigid>0)
    protein = protein->collapseRigidBonds(options.collapseRigid);

  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein->getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;

  //Initialize move
  Move* move;
  if(options.preventClashes){
    log("samplingStatus")<<"Using clash-avoiding move"<<endl;
    move = new ClashAvoidingMove( options.maxRotation,
                                  options.decreaseSteps,
                                  options.collisionCheck,
                                  options.projectConstraints );
  }else{
    log("samplingStatus")<<"Using nullspace move"<<endl;
    move = new NullspaceMove(RelativeTransitionOptions::getOptions()->maxRotation);

    if(options.decreaseSteps>0){
      log("samplingStatus")<<" .. with "<<options.decreaseSteps<<" decrease-steps"<<endl;
      move = new DecreaseStepMove(move, (unsigned int)options.decreaseSteps, options.decreaseFactor);
    }
  }
  move->setStepSize(options.stepSize);

  Selection resNetwork(options.residueNetwork);
  std::vector< std::tuple<Atom*, Atom*, double> > goal_distances =
      IO::readRelativeDistances(options.relativeDistances, protein);

  Direction* d1 = new LSNrelativeDirection(resNetwork, goal_distances);
  Direction* d2 = new RandomDirection(resNetwork,RelativeTransitionOptions::getOptions()->maxRotation);


  log() << "Total DOFs: " << protein->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs() << endl;fflush(stdout);
  log()<<"Number of rigid clusters: "<<protein->m_conf->m_numClusters;
  log()<<", biggest cluster: index "<<protein->m_conf->m_maxIndex<<" with "<<protein->m_conf->m_maxSize<<" atoms!"<<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidDihedrals() << " rigidified";
  log()<<" and " << ( protein->m_conf->getNullspace()->getNumDOFs()-
                      protein->m_conf->getNullspace()->getNumRigidDihedrals()) << " coordinated dihedrals" <<endl;
  log()<< protein->m_conf->getNullspace()->getNumRigidHBonds()<<" rigid out of "<<protein->getHBonds().size()<<" hydrogen bonds!"<<endl<<endl;


  log("samplingStatus")<<"Sampling ...\n";
  CTKTimer timer;
  timer.Reset();
  double start_time = timer.LastElapsedTime();
  gsl_vector* tmp1 = gsl_vector_alloc(protein->m_spanningTree->getNumDOFs());
  gsl_vector* tmp2 = gsl_vector_alloc(protein->m_spanningTree->getNumDOFs());
  gsl_vector* gradient = gsl_vector_alloc(protein->m_spanningTree->getNumDOFs());

  std::list<Configuration*> samples;
  samples.push_back(new Configuration(protein));
  double bias = options.biasToTarget;

  for(int i=0;i<options.samplesToGenerate;i++){
    double dist = dist_to_objective(goal_distances);
    if (dist<options.convergeDistance) break;

//    cout<<"Iteration "<<i<<endl;

    Configuration* seed = samples.back();
    d1->gradient(seed, seed, tmp1); //directed move
    d2->gradient(seed, seed, tmp2); //random move
    scale_gradient(tmp2, protein, options.maxRotation);
    double max_val1 = 0;
    double max_val2 = 0;
    for (int j=0;j<protein->m_spanningTree->getNumDOFs();j++){ //looking for the maximal rotation
      if ( fabs(gsl_vector_get(tmp1,j)) > max_val1)
        max_val1 = fabs(gsl_vector_get(tmp1,j));
      if ( fabs(gsl_vector_get(tmp2,j)) > max_val2)
        max_val2 = fabs(gsl_vector_get(tmp2,j));
    }
    if (max_val2==0 || bias==1){ //if the random gradient is null we don't use it
      for (int j = 0; j<protein->m_spanningTree->getNumDOFs();j++){
        gsl_vector_set(gradient, j, gsl_vector_get(tmp1,j));
      }
    }
    else{
      if (bias==0){ //only random
        for (int j = 0; j<protein->m_spanningTree->getNumDOFs();j++){
          gsl_vector_set(gradient, j, gsl_vector_get(tmp2,j));
        }
      }
      else { // we scale the random gradient so it has the 1/4 the max rotation as the directed gradient
        gsl_vector_scale(tmp2, max_val1 / (max_val2)*(1-bias)/bias);
        for (int j = 0; j<protein->m_spanningTree->getNumDOFs();j++){
          gsl_vector_set(gradient, j, gsl_vector_get(tmp1,j)+gsl_vector_get(tmp2,j));
        }
        //gsl_vector_scale_max_component(gradient, options.maxRotation);
      }
    }

    scale_gradient(gradient, protein, options.maxRotation);
    Configuration* new_conf = move->move(seed, gradient);
    samples.push_back(new_conf);

    log("samplingStatus")<<"> New structure: conf_"+std::to_string((long long)i)+".pdb"<<endl;
    string fname = "output/conf_"+std::to_string((long long)i)+".pdb";
    IO::writePdb(new_conf->updatedMolecule(), fname);
  }

  //Print final status
  double end_time = timer.ElapsedTime();
  log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(samples.size()-1)<<" valid samples\n";
  log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianTime<<" seconds\n";
  log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
  log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";
  log("samplingStatus")<< "Done"<<endl;
  delete d1;
  delete d2;


  return 0;
}


void scale_gradient(gsl_vector* gradient, Molecule* mol,double maxRotation)
{
  double max=1;
  for(int i=0;i<mol->m_spanningTree->getNumDOFs();i++){
    DOF* dof = mol->m_spanningTree->getDOF(i);
    int idx = dof->getIndex();
    double val = gsl_vector_get(gradient, idx);
    //double maxval = dof->getMaxPerturbation();
    if(fabs(val/maxRotation)>max)
      max = fabs(val/maxRotation);
  }
  gsl_vector_scale(gradient, 1/max);
}


double dist_to_objective(std::vector< std::tuple<Atom*, Atom*, double> > goal_distances)
{
  double d=0;
  for (int i=0; i< goal_distances.size();i++){
    Coordinate c1= get<0>(goal_distances[i])->m_position;
    Coordinate c2= get<1>(goal_distances[i])->m_position;

    d=d+fabs((sqrt((c1.x-c2.x)*(c1.x-c2.x)+(c1.y-c2.y)*(c1.y-c2.y)+(c1.z-c2.z)*(c1.z-c2.z))-get<2>(goal_distances[i])));
  }
  return d;
}
