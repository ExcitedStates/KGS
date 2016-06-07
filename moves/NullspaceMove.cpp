#include "NullspaceMove.h"

#include <math/gsl_helpers.h>
#include <gsl/gsl_vector_double.h>
#include "Logger.h"

using namespace std;

NullspaceMove::NullspaceMove(double maxRotation):
    m_maxRotation(maxRotation)
{

}


Configuration* NullspaceMove::performMove(Configuration* current, gsl_vector* gradient) {
  //current->computeCycleJacobianAndNullSpace();

  // Get atom positions at current
  //Molecule *protein = current->updatedMolecule();

  double currNorm = gsl_vector_length(gradient);
  log("dominik") << "Norm of gradient: " << currNorm << endl;

  // Project the gradient onto the null space of current
  gsl_vector *projected_gradient = gsl_vector_calloc(current->getNumDOFs());
  current->projectOnCycleNullSpace(gradient, projected_gradient);
  //protein->ProjectOnCycleNullSpace(gradient, projected_gradient);
  //gsl_vector *projected_gradient = projectOnNullspace( current, gradient, nullptr );

  double currProjNorm = gsl_vector_length(projected_gradient);
  log("dominik") << "Norm of projected gradient: " << currProjNorm << endl;

  Configuration *new_q = new Configuration(current);
  for (int i = 0; i < new_q->getNumDOFs(); ++i) {
    new_q->m_dofs[i] = current->m_dofs[i] + gsl_vector_get(projected_gradient, i);
    //new_q->m_dofs[i] = current->m_dofs[i] + min(m_stepSize, currNorm) * gsl_vector_get(projected_gradient, i);
  }

  gsl_vector_free(projected_gradient);
  return new_q;
}

gsl_vector* NullspaceMove::projectOnNullspace(Configuration* conf, gsl_vector* gradient, gsl_vector* ret)
{
  if(ret==nullptr) ret = gsl_vector_alloc(gradient->size);
  else assert(gradient->size==ret->size);

  // No cycles
  if(!conf->getNullspace()) {
    gsl_vector_memcpy(ret, gradient);
    return ret;
  }

  if( gradient->size >conf->getNullspace()->getNumDOFs() ) {
    // The input vectors contain all DOFs, however, the null space only contains DOFs in cycles.
    // Convert the DOFs in the input vectors to DOFs in cycles.

    //Fill reduced gradient
    gsl_vector *to_proj_short = gsl_vector_calloc(conf->getNullspace()->getNumDOFs());
    for (auto const& edge: conf->getMolecule()->m_spanning_tree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(to_proj_short,cycle_dof_id,gsl_vector_get(gradient, dof_id));
      }
    }

    //Project reduced gradient
    double normBefore = gsl_vector_length(to_proj_short);
    gsl_vector *after_proj_short = gsl_vector_calloc(conf->getNullspace()->getNumDOFs());
    conf->getNullspace()->ProjectOnNullSpace(to_proj_short, after_proj_short);
    double normAfter = gsl_vector_length(after_proj_short);

    //Scale so the norm is the same as before projection
    gsl_vector_scale(after_proj_short, normBefore/normAfter);

    //Convert back to full length DOFs vector
    for( auto const& edge: conf->getMolecule()->m_spanning_tree->Edges){
      int dof_id = edge->getDOF()->getIndex();
      int cycle_dof_id = edge->getDOF()->getCycleIndex();
      if ( cycle_dof_id!=-1 ) {
        gsl_vector_set(ret, dof_id,gsl_vector_get(after_proj_short,cycle_dof_id));
      }
      else if ( dof_id!=-1 ) {
        gsl_vector_set(ret, dof_id,gsl_vector_get(gradient,dof_id));
      }
    }
    gsl_vector_free(to_proj_short);
    gsl_vector_free(after_proj_short);
  } else {
    double normBefore = gsl_vector_length(gradient);
    conf->getNullspace()->ProjectOnNullSpace(gradient, ret);
    double normAfter = gsl_vector_length(ret);
    gsl_vector_scale(ret, normBefore/normAfter);
  }
  return ret;

}
