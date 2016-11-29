/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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

#include "FastClashAvoidingMove.h"

#include <math/gsl_helpers.h>
#include <math/SVDMKL.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <math/NullspaceSVD.h>
#include "Logger.h"

using namespace std;

FastClashAvoidingMove::FastClashAvoidingMove(
    double maxRotation,
    int trialSteps,
    const string& atomTypes,
    bool projectConstraints
) :
    m_maxRotation(maxRotation),
    m_trialSteps(trialSteps),
    m_collisionCheckAtomTypes(atomTypes),
    m_projectConstraints(projectConstraints)
//    m_maxRotation(SamplingOptions::getOptions()->maxRotation),
//    m_trialSteps(SamplingOptions::getOptions()->decreaseSteps),
//    m_collisionCheckAtomTypes(SamplingOptions::getOptions()->collisionCheck),
//    m_projectConstraints(SamplingOptions::getOptions()->projectConstraints)
{
  m_movesAccepted = 0;
  m_movesRejected = 0;
}

Configuration* FastClashAvoidingMove::performMove(Configuration* current, gsl_vector* gradient) {
//  enableLogger("clashBug");
  double currNorm = gsl_vector_length(gradient);
  log("dominik") << "Norm of gradient: " << currNorm << endl;

  // Project the gradient onto the null space of current
  gsl_vector *projected_gradient = gsl_vector_calloc(current->getNumDOFs());
  current->projectOnCycleNullSpace(gradient, projected_gradient);

  double currProjNorm = gsl_vector_length(projected_gradient);
  log("dominik") << "Norm of projected gradient: " << currProjNorm << endl;

  Configuration *new_q = new Configuration(current);
  for (int i = 0; i < new_q->getNumDOFs(); ++i)
    new_q->m_dofs[i] = formatRangeRadian( current->m_dofs[i] + gsl_vector_get(projected_gradient, i) );
  gsl_vector_free(projected_gradient);

  //If no clash return
  if (!new_q->updatedMolecule()->inCollision()) {
    m_movesAccepted++;
    return new_q;
  }
  //Else we use the clash avoiding/preventing move
  delete new_q; //added, otherwise memory leak with new definition below

  set< pair<Atom*,Atom*> > allCollisions, previousCollisions;

  //If resulting structure is in collision try for m_trialSteps times to avoid collision
  for (int trialStep = 0; trialStep < m_trialSteps; trialStep++) {

    //get all collisions at this configuration
    allCollisions = current->updatedMolecule()->getAllCollisions(m_collisionCheckAtomTypes);

    //Combine with collisions of previous trial
    for(auto const& prev_coll: previousCollisions){
      auto ait = allCollisions.find(prev_coll); //TODO: This shouldn't be necessary as we're using set. Verify
      if(ait==allCollisions.end()){
        allCollisions.insert(prev_coll);
      }
    }

    //This function overwrites necessary stuff in new_q
    Configuration* new_q = projectOnClashNullspace(current, gradient, allCollisions);

    // The new configuration is valid only if it is collision-free
    if (new_q->updatedMolecule()->inCollision()) {
      log("dominik") << "Rejected!" << endl;

      previousCollisions = allCollisions;
      if(trialStep==m_trialSteps-1){
        if(true) exit(-1);
        m_movesRejected++;
        return new_q;
      }

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

map<int,int> FastClashAvoidingMove::collectConstrainedDofMap(Configuration* conf, set< std::pair<Atom*,Atom*> >& allCollisions){
  //Associates general dof ids with constrained dof ids
  map<int,int> ret;

  //First add all cycle-DOFs
  for(auto const& edge: conf->getMolecule()->m_spanningTree->Edges){
    int cycle_dof_id = edge->getDOF()->getCycleIndex();
    int dof_id = edge->getDOF()->getIndex();
    if(cycle_dof_id>=0 && ret.count(dof_id)==0)
      ret[dof_id] = cycle_dof_id;
  }

  int nextidx = conf->getMolecule()->m_spanningTree->getNumCycleDOFs();

  //Add all clash-DOFs
  for(auto const& coll: allCollisions) {
    KinVertex *vertex1 = coll.first->getRigidbody()->getVertex();
    KinVertex *vertex2 = coll.second->getRigidbody()->getVertex();
    KinVertex *common_ancestor = conf->getMolecule()->m_spanningTree->findCommonAncestor(vertex1, vertex2);

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

Configuration* FastClashAvoidingMove::projectOnClashNullspace(
    Configuration *conf,
    gsl_vector *gradient,
    set<std::pair<Atom *, Atom *> > &collisions
){
//  log("clashBug")<<"projectOnClashNullspace(..)"<<endl;

  //Associates general dof ids with constrained dof ids.
  map<int,int> constrainedDofMap = collectConstrainedDofMap(conf, collisions);

  //Transfer the parts of the gradient that are in clash or constraint cycles.
  gsl_vector* reducedGradient = gsl_vector_alloc(constrainedDofMap.size());
  for(auto const& general_clash_pair: constrainedDofMap){
    double generalValue = gsl_vector_get(gradient, general_clash_pair.first);
//    log("clashBug")<<"> dof "<<general_clash_pair.first<<" (general) / "<<general_clash_pair.second<<" (constrained) = "<<generalValue<<endl;
    gsl_vector_set(reducedGradient, general_clash_pair.second, generalValue);
  }

  //Compute clash-avoiding jacobian, svd, and nullspace
  //Todo: This could be optimized, we only have to compute the Jacobian the first time...
  gsl_matrix* clashJac = computeClashAvoidingJacobian(conf, constrainedDofMap, collisions);
  SVD* clashSVD = SVD::createSVD(clashJac);//new SVDMKL(clashAvoidingJacobian);
  Nullspace* clashNullSpace = new NullspaceSVD(clashSVD);
  clashNullSpace->updateFromMatrix();

  //Project reducedGradient
  double normBefore = gsl_vector_length(reducedGradient);
  clashNullSpace->projectOnNullSpace(reducedGradient, reducedGradient);
  double normAfter = gsl_vector_length(reducedGradient);
//  log("clashBug")<<"> normBefore: "<<normBefore<<endl;
//  log("clashBug")<<"> normAfter:  "<<normAfter<<endl;

  //Scale so the length matches the one before
  if(normAfter>0.00000001)
    gsl_vector_scale(reducedGradient, normBefore/normAfter);

  //Transfer to general dofs again
  gsl_vector* projected_gradient = gsl_vector_copy(gradient);
  for(auto const& general_clash_pair: constrainedDofMap){
    double constrainedValue = gsl_vector_get(reducedGradient, general_clash_pair.second);
    gsl_vector_set(projected_gradient, general_clash_pair.first, constrainedValue);
  }

  //Clean up
  gsl_vector_free(reducedGradient);
  gsl_matrix_free(clashJac);

  double currProjNorm = gsl_vector_length(projected_gradient);

  Configuration* new_q = new Configuration(conf);
//  log("dominik")<<"Clash trial "<<trialStep<<", Norm of projected gradient: "<<currProjNorm<<endl;
  for (int i = 0; i < new_q->getNumDOFs(); ++i)
    new_q->m_dofs[i] = formatRangeRadian( conf->m_dofs[i] + gsl_vector_get(projected_gradient, i) );
  gsl_vector_free(projected_gradient);

  new_q->m_usedClashPrevention = true;
  new_q->m_clashFreeDofs = new_q->getNumDOFs() - constrainedDofMap.size() + clashNullSpace->getNullspaceSize();

  //Clean up
  delete clashNullSpace;
  delete clashSVD;

  return new_q;
}

gsl_matrix* FastClashAvoidingMove::computeClashAvoidingJacobian(
    Configuration* conf,
    map<int,int>& dofMap,
    set< std::pair<Atom*,Atom*> >& collisions
) {
  //The clash Jacobian is the regular Jacobian's constraints, plus one constraint per pair of clashing atoms

  //Clashes can occur also for previously free dihedrals!
  //Therefore, we use the full set of dihedrals to determine this matrix!
  conf->updateMolecule();
  //Check that the correct cycle jacobian is used.
  gsl_matrix *cycleJac = conf->getCycleJacobian();
  int numCollisions = collisions.size();
  int rowNum = cycleJac->size1 + numCollisions;
  //int colNum = conf->getMolecule()->m_spanningTree->getNumDOFs();

  //No longer uses all dihedrals as it messes with scaling
  int colNum = dofMap.size();

  gsl_matrix *ret = gsl_matrix_calloc(rowNum, colNum);

  //Copy cycle-jacobian into ret
  for (int r = 0; r < cycleJac->size1; r++) {
    for (int c = 0; c < cycleJac->size2; c++) {
      gsl_matrix_set(ret, r, c, gsl_matrix_get(cycleJac, r, c));
    }
  }

//  //Convert the cycle Jacobian to a full Jacobian
//  //Columns correspond to cycle_dof_ids
//  if(m_projectConstraints){
//    map<unsigned int, KinVertex*>::iterator vit;
//
//    for(auto const& edge: conf->getMolecule()->m_spanningTree->Edges){
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

    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2

    Math3D::Vector3 clashNormal = p2-p1;
    clashNormal.getNormalized(clashNormal);

    //Vertices
    KinVertex* vertex1 = atom1->getRigidbody()->getVertex();
    KinVertex* vertex2 = atom2->getRigidbody()->getVertex();
    KinVertex* common_ancestor = conf->getMolecule()->m_spanningTree->findCommonAncestor(vertex1, vertex2);

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
