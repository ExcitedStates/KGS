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
// Created by Dominik Budday on 31.07.18.
//

#include "LSNclashAvoidingMove.h"

#include <math/gsl_helpers.h>
#include <math/SVDMKL.h>
#include <math/NullspaceSVD.h>
#include "Logger.h"
#include "directions/LSNrelativeDirection.h"

using namespace std;

LSNclashAvoidingMove::LSNclashAvoidingMove(LSNrelativeDirection *direction,
    double maxRotation,
    int trialSteps,
    const std::string& atomTypes,
    bool projectConstraints
) :
    Move(maxRotation),
    m_direction(direction),
    m_trialSteps(trialSteps),
    m_collisionCheckAtomTypes(atomTypes),
    m_projectConstraints(projectConstraints)
{
  m_movesAccepted = 0;
  m_movesRejected = 0;
}

Configuration* LSNclashAvoidingMove::performMove(Configuration* current, gsl_vector* projected_gradient) {
//  current->computeJacobians();

  // Get atom positions at current
  Molecule *protein = current->updatedMolecule();

  //clash prevention technique
  bool usedClashJacobian = false;
  set< pair<Atom*,Atom*> > allCollisions, previousCollisions;

  // Create new configuration
  Configuration* new_q = new Configuration(current);
  new_q->m_clashFreeDofs = new_q->getNumDOFs() - new_q->getMolecule()->m_spanningTree->getNumCycleDOFs() +
                           new_q->getNullspace()->getNullspaceSize();

  //If resulting structure is in collision try scaling down the gradient
  for (int trialStep = 0; trialStep <= m_trialSteps; trialStep++) {
    //Scale max entry
    if(m_scale)
      gsl_vector_scale_max_component(projected_gradient,m_maxRotation);
    double currProjNorm = gsl_vector_length(projected_gradient);

    for (int i = 0; i < new_q->getNumDOFs(); ++i) {
      new_q->m_dofs[i] = formatRangeRadian( current->m_dofs[i] + gsl_vector_get(projected_gradient, i) ); //added as transformation always starts from original coordinates
    }

    new_q->getMolecule()->forceUpdateConfiguration(new_q);

    // The new configuration is valid only if it is collision-free
    if (new_q->updatedMolecule()->inCollision()) {
      cout<<"Rejected, computing clash-free move"<<endl;
      log("planner") << "Rejected!" << endl;
      m_movesRejected++;

//      log("planner")<<"Now computing a clash free m_direction!"<<endl;
      allCollisions = protein->getAllCollisions(m_collisionCheckAtomTypes);//get all collisions at this configuration

      for(auto const& prev_coll: previousCollisions){//combine collisions
        auto ait = allCollisions.find(prev_coll);
        if(ait==allCollisions.end()){
          allCollisions.insert(prev_coll);
        }
      }

      //Now we have the set of collisions we use as additional constraints in the new Jacobian
      current->updateMolecule();

      if(!usedClashJacobian) {//first clash-avoiding trial, recompute Jacobian matrices
        current->getCycleJacobian();
        usedClashJacobian = true; //flag to recompute Jacobian
      }

      gsl_matrix* clashAvoidingJacobian = computeClashAvoidingJacobian( current, allCollisions);
      SVD* clashAvoidingSVD  = new SVDMKL(clashAvoidingJacobian);
      Nullspace* clashAvoidingNullSpace = new NullspaceSVD(clashAvoidingSVD);
      clashAvoidingNullSpace->updateFromMatrix();
      ///Computing a new Least-Square solution with the new clash-free nullspace
      m_direction->computeGradientExternalN(clashAvoidingNullSpace->getBasis(), current, nullptr, projected_gradient);

      //Scale max entry
      if(m_scale)
        gsl_vector_scale_max_component(projected_gradient,m_maxRotation);

      new_q->m_usedClashPrevention = true;
      new_q->m_clashFreeDofs = clashAvoidingNullSpace->getNullspaceSize();

      delete clashAvoidingNullSpace;
      delete clashAvoidingSVD;
      gsl_matrix_free(clashAvoidingJacobian);

      previousCollisions=allCollisions; //save collisions in case the new one is rejected again

    } else {//collision free
      m_movesAccepted++;
      return new_q;
    }
  }//end steps

  return new_q;
}

gsl_matrix* LSNclashAvoidingMove::computeClashAvoidingJacobian(Configuration* conf, set< std::pair<Atom*,Atom*> >& allCollisions) {

  //The clash Jacobian is the regular Jacobian's constraints, plus one constraint per pair of clashing atoms
  int numCollisions = allCollisions.size();

  //Clashes can occur also for previously free dihedrals!
  //Therefore, we use the full set of dihedrals to determine this matrix!

  int rowNum = conf->getCycleJacobian()->size1 + numCollisions;
  int colNum = conf->getMolecule()->m_spanningTree->getNumDOFs();

  gsl_matrix *ret = gsl_matrix_calloc(rowNum, colNum);

  //Convert the cycle Jacobian to a full Jacobian
  //Columns correspond to cycle_dof_ids
  if (m_projectConstraints) {

    gsl_matrix *cycleJac = conf->getCycleJacobian();

    for (auto const &edge: conf->getMolecule()->m_spanningTree->m_edges) {
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if (cycle_dof_id != -1) {
        for (int i = 0; i != conf->getCycleJacobian()->size1; i++) {
          gsl_matrix_set(ret, i, dof_id, gsl_matrix_get(cycleJac, i, cycle_dof_id));
        }
      }
    }
  } // else: we maintain a zero-valued Jacobian to not consider the constraints (only for testing of constraint limitations)

  int i = conf->getCycleJacobian()->size1;//start entries after all cycles

  //Quick trick for plotting the clash-cycle network
//	int consCounter = 1;
//	vector<Atom*>::iterator ait;
//	for(ait=protein->Atom_list.begin(); ait != protein->Atom_list.end(); ait++){
//		(*ait)->m_assignedBiggerRB_id = consCounter;
//	}

  for (auto const &coll: allCollisions) {
    Atom *atom1 = coll.first;
    Atom *atom2 = coll.second;
    log("planner") << "Using clash constraint for atoms: " << atom1->getId() << " " << atom2->getId() << endl;

    Coordinate p1 = atom1->m_position; //end-effector, position 1
    Coordinate p2 = atom2->m_position; //end-effector, position 2

    Math3D::Vector3 clashNormal = p2 - p1;
    clashNormal.getNormalized(clashNormal);

    //Vertices
    KinVertex *vertex1 = atom1->getRigidbody()->getVertex();
    KinVertex *vertex2 = atom2->getRigidbody()->getVertex();
    KinVertex *common_ancestor = conf->getMolecule()->m_spanningTree->findCommonAncestor(vertex1, vertex2);

    // trace back until the common ancestor from vertex1
    while (vertex1 != common_ancestor) {
      KinVertex *parent = vertex1->m_parent;
      KinEdge *p_edge = parent->findEdge(vertex1);

      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id != -1) { // this edge is a DOF

//        Math3D::Vector3 derivativeP1 = ComputeJacobianEntry(p_edge->getBond()->Atom1->m_position,p_edge->getBond()->m_atom2->m_position,p1);
        Math3D::Vector3 derivativeP1 = p_edge->getDOF()->getDerivative(p1);
        double jacobianEntryClash = dot(clashNormal, derivativeP1);

        gsl_matrix_set(ret, i, dof_id, jacobianEntryClash); //set: Matrix, row, column, what to set
//				log("planner")<<"Setting clash entry on left branch at "<<dof_id<<endl;
      }

      //Quick trick for plotting the clash cycles (don't use in real sampling)
//			for (ait=vertex1->Rb_ptr->Atoms.begin(); ait !=vertex1->Rb_ptr->Atoms.end(); ait++){
//				Atom* atom = (*ait);
//				atom->m_assignedBiggerRB_id = consCounter;
//			}
      vertex1 = parent;
    }

    // trace back until the common ancestor from vertex2
    while (vertex2 != common_ancestor) {
      KinVertex *parent = vertex2->m_parent;
      KinEdge *p_edge = parent->findEdge(vertex2);

      int dof_id = p_edge->getDOF()->getIndex();
      if (dof_id != -1) { // this edge is a DOF

//        Math3D::Vector3 derivativeP2 = ComputeJacobianEntry(
//            p_edge->getBond()->Atom1->m_position,
//            p_edge->getBond()->m_atom2->m_position,
//            p2); //b

        Math3D::Vector3 derivativeP2 = p_edge->getDOF()->getDerivative(p2);
        double jacobianEntryClash = -dot(clashNormal, derivativeP2);

        gsl_matrix_set(ret, i, dof_id, jacobianEntryClash); //set: Matrix, row, column, what to set
      }

      vertex2 = parent;
    }
    ++i;
  }

  return ret;
}