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

#include <moves/NullspaceMove.h>
#include <moves/ClashAvoidingMove.h>
#include <metrics/Dihedral.h>
#include <directions/RandomDirection.h>
#include <directions/DihedralDirection.h>
#include <directions/MSDDirection.h>
#include <directions/LSNullspaceDirection.h>
#include <directions/BlendedDirection.h>
#include <moves/DecreaseStepMove.h>
#include <metrics/RMSDnosuper.h>

using namespace std;

extern double jacobianTime;
extern double rigidityTime;
extern double selectNodeTime;

void scale_gradient(gsl_vector* gradient, Molecule* mol);

int main( int argc, char* argv[] ) {
  enableLogger("default");
  enableLogger("samplingStatus");

  //SamplingOptions options(argc,argv);
  SamplingOptions::createOptions(argc, argv);

  SamplingOptions &options = *(SamplingOptions::getOptions());

  if (loggerEnabled("samplingStatus")) {
    enableLogger("so");//SamplingOptions
    options.print();
  }

  // Set seed
  srand(options.seed);

  Selection movingResidues(options.residueNetwork);
  Molecule* protein = IO::readPdb(
      options.initialStructureFile,
      movingResidues,
      options.extraCovBonds,
      options.roots,
      options.hydrogenbondMethod,
      options.hydrogenbondFile
  );
  protein->setCollisionFactor(options.collisionFactor);

//  if(options.hydrogenbondMethod=="user")
//    IO::readHbonds( &protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="rnaview")
//    IO::readHbonds_rnaview( &protein, options.hydrogenbondFile, options.annotationFile.empty() );
//  else if(options.hydrogenbondMethod=="first" || options.hydrogenbondMethod=="FIRST")
//    IO::readHbonds_first( &protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="vadar")
//    IO::readHbonds_vadar( &protein, options.hydrogenbondFile );
//  else if(options.hydrogenbondMethod=="dssr")
//    IO::readHbonds_dssr( &protein, options.hydrogenbondFile );


  log("samplingStatus")<<"Molecule has:"<<endl;
  log("samplingStatus")<<"> "<<protein->getAtoms().size() << " atoms" << endl;
  log("samplingStatus")<<"> "<<protein->getInitialCollisions().size()<<" initial collisions"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->m_cycleAnchorEdges.size()<<" hydrogen bonds"<<endl;
  log("samplingStatus")<<"> "<<protein->m_spanningTree->getNumDOFs() << " DOFs of which " << protein->m_spanningTree->getNumCycleDOFs() << " are cycle-DOFs\n" << endl;


  //Initialize metric
  metrics::Metric* metric = nullptr;
  try {
    Selection metricSelection(options.metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="rmsd") 		    metric = new metrics::RMSD(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="rmsdnosuper") metric = new metrics::RMSDnosuper(metricSelection);
    if(SamplingOptions::getOptions()->metric_string=="dihedral")    metric = new metrics::Dihedral(metricSelection);
  }catch(std::runtime_error& error) {
    cerr<<error.what()<<endl;
    exit(-1);
  }

  //Initialize move
  Move* move;
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
  Direction* direction;
  bool blendedDir = false;
  if(options.gradient == 0)
    direction = new RandomDirection(movingResidues);
  else if(options.gradient == 1)
    direction = new DihedralDirection(movingResidues);
  else if(options.gradient == 2){
    BlendedDirection* m_direction = new BlendedDirection();
    m_direction->addDirection(new DihedralDirection(movingResidues),0);
    m_direction->addDirection(new RandomDirection(movingResidues,SamplingOptions::getOptions()->maxRotation), 1);
    direction = m_direction;
    blendedDir = true;
  }
  else if(options.gradient == 3)
    direction = new MSDDirection(movingResidues);
  else if(options.gradient == 4){
    BlendedDirection* m_direction = new BlendedDirection();
    m_direction->addDirection(new MSDDirection(movingResidues),0);
    m_direction->addDirection(new RandomDirection(movingResidues,SamplingOptions::getOptions()->maxRotation), 1);
    direction = m_direction;
    blendedDir = true;
  }
  else if(options.gradient <= 5)
    direction = new LSNullspaceDirection(movingResidues);


  if(options.saveData > 0){

    log() << "Total DOFs: " << protein->m_spanningTree->getNumDOFs() << ", Cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs() << endl;fflush(stdout);

    if(options.saveData > 1){
      string out = options.workingDirectory + "output/" + protein->getName() + "_q_0.txt";
      IO::writeQ(protein, protein->m_conf, out);
    }

    log()<<"Number of rigid clusters: "<<protein->m_conf->m_numClusters;
    log()<<", biggest cluster: index "<<protein->m_conf->m_maxIndex<<" with "<<protein->m_conf->m_maxSize<<" atoms!"<<endl;
    log()<< protein->m_conf->getNullspace()->getNumRigidDihedrals() << " rigidified";
    log()<<" and " << ( protein->m_conf->getNullspace()->getNumDOFs()-
                        protein->m_conf->getNullspace()->getNumRigidDihedrals()) << " coordinated dihedrals" <<endl;
    log()<< protein->m_conf->getNullspace()->getNumRigidHBonds()<<" rigid out of "<<protein->getHBonds().size()<<" hydrogen bonds!"<<endl<<endl;

    double bigRad = options.stepSize*4/3;
    double lilRad = bigRad/2;

    log("samplingStatus")<<"Sampling ...\n";
    CTKTimer timer;
    timer.Reset();
    double start_time = timer.LastElapsedTime();

    int clashes = 0;

    gsl_vector* gradient = gsl_vector_alloc(protein->m_spanningTree->getNumDOFs());
    std::list<Configuration*> samples;
    samples.push_back(new Configuration(protein));
    for(int i=0;i<options.samplesToGenerate;i++){
      Configuration* seed = samples.back();
      direction->gradient(seed, nullptr, gradient);
      //gsl_vector_scale_max_component(gradient, options.maxRotation);

      // Scale gradient so move is in Poisson disc
      Configuration *new_conf = move->move(seed, gradient); //Perform move
      double dist = metric->distance(new_conf, seed);
      int scaleAttempts = 0;
      while( dist<lilRad || dist>bigRad){
        if(++scaleAttempts==5) break;
        double gradientScale = (bigRad+lilRad)/(2.0*dist);
        gsl_vector_scale(gradient, gradientScale);
        delete new_conf;
        new_conf = move->move(seed, gradient);
        dist = metric->distance(new_conf, seed);
      }

//      scale_gradient(gradient, &protein);
//      gsl_vector_scale(gradient, options.stepSize);
//      Configuration* new_conf = move->move(seed, gradient);

      if(new_conf->updatedMolecule()->inCollision()) {
        clashes++;
        i--;
      } else {
        log("samplingStatus")<<"Accepted conformation "<<i<<endl;
        IO::writePdb(new_conf->updatedMolecule(), "output/conf_" + std::to_string((long long) i) + ".pdb");
        samples.push_back(new_conf);
      }
    }


    //Print final status
    double end_time = timer.ElapsedTime();
    log("samplingStatus")<< "Took "<<(end_time-start_time)<<" seconds to generate "<<(samples.size()-1)<<" valid samples\n";
    log("samplingStatus")<< "Jacobian and null space computation took "<<jacobianTime<<" seconds\n";
    log("samplingStatus")<< "Rigidity analysis took "<<rigidityTime<<" seconds\n";
    log("samplingStatus")<< "Node selection took "<<selectNodeTime<<" seconds\n";
    log("samplingStatus")<< "Clashes: "<<clashes<<endl;
//    log("samplingStatus")<< planner->initialRebuildsAccepted<<" initial rebuild perturbations accepted (RI_ACC)"<<endl;
//    log("samplingStatus")<< planner->initialRebuildsRejected<<" initial rebuild perturbations rejected for clashing (RI_CLASH)"<<endl;
//    log("samplingStatus")<< planner->rebuildsAccepted<<" rebuild perturbations accepted (R_ACC)"<<endl;
//    log("samplingStatus")<< planner->rebuildsRejected<<" rebuild perturbations rejected for clashing (R_CLASH)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sAccepted<<" nullspace perturbations accepted (NS_ACC)"<<endl;
//    log("samplingStatus")<< planner->getNullspace()sRejected<<" nullspace perturbations rejected for clashing (NS_CLASH)"<<endl;
    log("samplingStatus")<<"MCMC-planner___: "<<samples.size()<<" open samples and 0";
    log("samplingStatus")<<" closed samples on termination ";
    log("samplingStatus")<<"("<<samples.size()<<" total)"<<endl;
    log("samplingStatus")<<"MCMC-planner___: Rejects from clash:          "<<clashes<<endl;
    log("samplingStatus")<<"MCMC-planner___: Rejects from tree-collision: 0"<<endl;

    if(options.saveData > 0){
      log("samplingStatus")<<"Creating trajectory"<<endl;
    }
  }
  log("samplingStatus")<<"Done"<<endl;
  //Clean up
  delete direction;


  return 0;
}


void scale_gradient(gsl_vector* gradient, Molecule* mol)
{
  double factor = 1.0;

  for(int i=0;i<mol->m_spanningTree->getNumDOFs();i++){
    DOF* dof = mol->m_spanningTree->getDOF(i);
    int idx = dof->getIndex();
    double val = gsl_vector_get(gradient, idx);
    double maxval = dof->getMaxPerturbation();
    if(fabs(maxval/val)<factor)
      factor = fabs(maxval/val);
  }
  gsl_vector_scale(gradient, factor);
}

