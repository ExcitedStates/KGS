//
// Created by Dominik Budday on 15.03.16.
//

#include "ClashAvoidingMove.h"

#include <math/gsl_helpers.h>
#include <math/MKLSVD.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include "Logger.h"

using namespace std;

ClashAvoidingMove::ClashAvoidingMove() :
    m_maxRotation(SamplingOptions::getOptions()->maxRotation),
    m_trialSteps(SamplingOptions::getOptions()->decreaseSteps),
    m_collisionCheckAtomTypes(SamplingOptions::getOptions()->collisionCheck),
    m_projectConstraints(SamplingOptions::getOptions()->projectConstraints)
{
  m_movesAccepted = 0;
  m_movesRejected = 0;
}

Configuration* ClashAvoidingMove::performMove(Configuration* current, gsl_vector* gradient) {
//  current->computeJacobians();
  //cout<<"ClashAvoidingMove::performMove(..)"<<endl;
  //cout<<"gradient: ";
  //for(int i=0;i<10;i++)
  //  cout<<gsl_vector_get(gradient, i)<<" ";
  //cout<<endl;
  //gsl_vector_cout(gradient);

  // Get atom positions at current
  Molecule *protein = current->updatedMolecule();

  double currNorm = gsl_vector_length(gradient);
  double targetNorm = currNorm;//*m_stepSize;

  log("dominik") << "Norm of gradient: " << currNorm << endl;
  //cout << "Norm of gradient: " << currNorm << endl;

  // Project the gradient onto the null space of current
  gsl_vector *projected_gradient = gsl_vector_calloc(protein->totalDofNum());
  protein->ProjectOnCycleNullSpace(gradient, projected_gradient);

  Configuration* new_q = new Configuration(current);
  std::copy( projected_gradient->data, projected_gradient->data+projected_gradient->size, new_q->m_dofs );
  gsl_vector_free(projected_gradient);

  //If no clash return
  if(!new_q->updatedMolecule()->inCollision())
    return new_q;

  set< pair<Atom*,Atom*> > allCollisions, previousCollisions;
  //int clashAvoidingDOFs = -1;

  //If resulting structure is in collision try scaling down the gradient
  for (int trialStep = 0; trialStep < m_trialSteps; trialStep++) {
    //cout<<"ClashAvoidingMove - trial: "<<trialStep<<" / "<<m_trialSteps<<endl;

    //get all collisions at this configuration
    allCollisions = protein->getAllCollisions(m_collisionCheckAtomTypes);

    //Combine with collisions of previous trial
    for(auto const& prev_coll: previousCollisions){
      auto ait = allCollisions.find(prev_coll); //TODO: This shouldn't be necessary as we're using set. Verify
      if(ait==allCollisions.end()){
        allCollisions.insert(prev_coll);
      }
    }

    projected_gradient = projectOnClashNullspace(current, gradient, allCollisions);
    double currProjNorm = gsl_vector_length(projected_gradient);
    log("dominik")<<"Norm of projected gradient: "<<currProjNorm<<endl;
    new_q = new Configuration(current);
    std::copy(
        projected_gradient->data,
        projected_gradient->data+new_q->getNumDOFs(),
        new_q->m_dofs );
    gsl_vector_free(projected_gradient);


    // The new configuration is valid only if it is collision-free
    if (new_q->updatedMolecule()->inCollision()) {
      log("dominik") << "Rejected!" << endl;
      m_movesRejected++;

      previousCollisions = allCollisions;
      if(trialStep==m_trialSteps-1) return new_q;

      delete new_q;

    } else {//collision free
      log("dominik") << "Accepted!" << endl;
      m_movesAccepted++;

      return new_q;
    }
  }//end steps

  cerr<<"ClashAvoidingMove::performMove - should not reach this point"<<endl;
  throw "ClashAvoidingMove::performMove - should not reach this point";
}

map<int,int> ClashAvoidingMove::collectConstrainedDofMap(Configuration* conf, set< std::pair<Atom*,Atom*> >& allCollisions){
  //Associates general dof ids with constrained dof ids
  map<int,int> ret;

  //First add all cycle-DOFs
  for(auto const& edge: conf->getMolecule()->m_spanning_tree->Edges){
    int cycle_dof_id = edge->getDOF()->getCycleIndex();
    int dof_id = edge->getDOF()->getIndex();
    if(cycle_dof_id>=0 && ret.count(dof_id)==0)
      ret[dof_id] = cycle_dof_id;
  }

  int nextidx = conf->getMolecule()->m_spanning_tree->getNumCycleDOFs();

  //Add all clash-DOFs
  for(auto const& coll: allCollisions) {
    KinVertex *vertex1 = coll.first->getRigidbody()->getVertex();
    KinVertex *vertex2 = coll.second->getRigidbody()->getVertex();
    KinVertex *common_ancestor = conf->getMolecule()->m_spanning_tree->findCommonAncestor(vertex1, vertex2);

    // trace back until the common ancestor from vertex1
    while (vertex1 != common_ancestor) {
      KinEdge* p_edge = vertex1->m_parent->findEdge(vertex1);

      int dof_id = p_edge->getDOF()->getIndex();
      if(ret.count(dof_id)==0)
        ret[dof_id] = nextidx++;

      vertex1 = vertex1->m_parent;
    }

    // trace back until the common ancestor from vertex2
    while (vertex2 != common_ancestor) {
      KinEdge* p_edge = vertex2->m_parent->findEdge(vertex2);

      int dof_id = p_edge->getDOF()->getIndex();
      if(ret.count(dof_id)==0)
        ret[dof_id] = nextidx++;

      vertex2 = vertex2->m_parent;
    }
  }
  return std::move(ret);
}

gsl_vector* ClashAvoidingMove::projectOnClashNullspace(
    Configuration *conf,
    gsl_vector *gradient,
    set<std::pair<Atom *, Atom *> > &collisions
){

  //Associates general dof ids with constrained dof ids.
  map<int,int> constrainedDofMap = collectConstrainedDofMap(conf, collisions);

  //Transfer the parts of the gradient that are in clash or constraint cycles.
  gsl_vector* reducedGradient = gsl_vector_alloc(constrainedDofMap.size());
  for(auto const& general_clash_pair: constrainedDofMap){
    double generalValue = gsl_vector_get(gradient, general_clash_pair.first);
    gsl_vector_set(reducedGradient, general_clash_pair.second, generalValue);
  }

  //Compute clash-avoiding jacobian, svd, and nullspace
  gsl_matrix* clashJac = computeClashAvoidingJacobian(conf, constrainedDofMap, collisions);
  SVD* clashSVD = SVD::createSVD(clashJac);//new MKLSVD(clashAvoidingJacobian);
  Nullspace* clashNullSpace = new Nullspace(clashSVD);
  clashNullSpace->UpdateFromMatrix();

  //Project reducedGradient
  double normBefore = gsl_vector_length(reducedGradient);
  clashNullSpace->ProjectOnNullSpace(reducedGradient, reducedGradient);
  double normAfter = gsl_vector_length(reducedGradient);

  //Scale so the length matches the one before
  gsl_vector_scale(reducedGradient, normBefore/normAfter);

  //Transfer to general dofs again
  gsl_vector* ret = gsl_vector_copy(gradient);
  for(auto const& general_clash_pair: constrainedDofMap){
    double constrainedValue = gsl_vector_get(reducedGradient, general_clash_pair.second);
    gsl_vector_set(ret, general_clash_pair.first, constrainedValue);
  }

  //Clean up
  gsl_vector_free(reducedGradient);
  gsl_matrix_free(clashJac);
  delete clashNullSpace;
  delete clashSVD;

  return ret;
}

gsl_matrix* ClashAvoidingMove::computeClashAvoidingJacobian(
    Configuration* conf,
    map<int,int>& dofMap,
    set< std::pair<Atom*,Atom*> >& collisions
)
{
  //The clash Jacobian is the regular Jacobian's constraints, plus one constraint per pair of clashing atoms

  //Clashes can occur also for previously free dihedrals!
  //Therefore, we use the full set of dihedrals to determine this matrix!

  gsl_matrix* cycleJac = conf->getCycleJacobian();
  int numCollisions = collisions.size();
  int rowNum = cycleJac->size1 + numCollisions;
  //int colNum = conf->getMolecule()->m_spanning_tree->getNumDOFs();

  //No longer uses all dihedrals as it messes with scaling
  int colNum = dofMap.size();

  gsl_matrix* ret = gsl_matrix_calloc(rowNum, colNum);

  //Copy cycle-jacobian into ret
  for(int r=0;r<cycleJac->size1;r++){
    for(int c=0;c<cycleJac->size2;c++){
      gsl_matrix_set(ret, r, c, gsl_matrix_get(cycleJac,r,c));
    }
  }
//  //Convert the cycle Jacobian to a full Jacobian
//  //Columns correspond to cycle_dof_ids
//  if(m_projectConstraints){
//    map<unsigned int, KinVertex*>::iterator vit;
//
//    for(auto const& edge: conf->getMolecule()->m_spanning_tree->Edges){
//      int dof_id = edge->getDOF()->getIndex();
//      int cycle_dof_id = edge->getDOF()->getCycleIndex();
//      if ( cycle_dof_id!=-1 ) {
//        for( int i=0; i!=conf->getCycleJacobian()->size1; i++){
//          gsl_matrix_set(ret, i, cycle_dof_id, gsl_matrix_get(conf->getCycleJacobian(),i,cycle_dof_id));
//        }
//      }
//    }
//  } // else: we maintain a zero-valued Jacobian to not consider the constraints (only for testing of constraint limitations)

  int r=cycleJac->size1;//start entries after all cycles
  int idx = cycleJac->size2;

  for(auto const& coll: collisions){
    Atom* atom1 = coll.first;
    Atom* atom2 = coll.second;
    log("dominik") << "Using clash constraint for atoms: "<<atom1->getId() << " " << atom2->getId() << endl;

    Coordinate p1 = atom1->m_Position; //end-effector, position 1
    Coordinate p2 = atom2->m_Position; //end-effector, position 2

    Math3D::Vector3 clashNormal = p2-p1;
    clashNormal.getNormalized(clashNormal);

    //Vertices
    KinVertex* vertex1 = atom1->getRigidbody()->getVertex();
    KinVertex* vertex2 = atom2->getRigidbody()->getVertex();
    KinVertex* common_ancestor = conf->getMolecule()->m_spanning_tree->findCommonAncestor(vertex1, vertex2);

    // trace back until the common ancestor from vertex1
    while ( vertex1 != common_ancestor ) {
      KinVertex* parent = vertex1->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex1);

      //Locate constrained dof id
      int dof_id = p_edge->getDOF()->getIndex();
      int constrained_dof_id = dofMap[dof_id];

      Math3D::Vector3 derivativeP1 = p_edge->getDOF()->getDerivative(p1);
      double jacobianEntryClash = dot(clashNormal, derivativeP1);

      gsl_matrix_set(ret,r,constrained_dof_id,jacobianEntryClash); //set: Matrix, row, column, what to set

      vertex1 = parent;
    }

    // trace back until the common ancestor from vertex2
    while ( vertex2 != common_ancestor ) {
      KinVertex* parent = vertex2->m_parent;
      KinEdge* p_edge = parent->findEdge(vertex2);

      int dof_id = p_edge->getDOF()->getIndex();
      int constrained_dof_id = dofMap[dof_id];

      Math3D::Vector3 derivativeP2 = p_edge->getDOF()->getDerivative(p2);

      double jacobianEntryClash = - dot(clashNormal, derivativeP2);

      gsl_matrix_set(ret,r,constrained_dof_id,jacobianEntryClash); //set: Matrix, row, column, what to set

      vertex2 = parent;
    }

    r++; //Next collision
  }

  return ret;
}
