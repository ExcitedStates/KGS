#include "NullspaceMove.h"

#include <math/gsl_helpers.h>
#include "Logger.h"

using namespace std;

NullspaceMove::NullspaceMove(double maxRotation):
    m_maxRotation(maxRotation)
{

}

Configuration* NullspaceMove::performMove(Configuration* current, gsl_vector* gradient) {
  current->computeCycleJacobianAndNullSpace();

  // Get atom positions at current
  Molecule *protein = current->updatedMolecule();

  double currNorm = gsl_vector_length(gradient);
  log("dominik") << "Norm of gradient: " << currNorm << endl;

  // Project the gradient onto the null space of current
  gsl_vector *projected_gradient = gsl_vector_calloc(protein->totalDofNum());
  protein->ProjectOnCycleNullSpace(gradient, projected_gradient);

  double currProjNorm = gsl_vector_length(projected_gradient);
  log("dominik") << "Norm of projected gradient: " << currProjNorm << endl;

  // Normalize projected_gradient
  if (currProjNorm > 0.001)
    gsl_vector_normalize(projected_gradient);

  gsl_vector_scale_max_component(projected_gradient, m_maxRotation);

  Configuration *new_q = new Configuration(current);

  for (int i = 0; i < new_q->getNumDOFs(); ++i) {
    //TODO: Optimize using gsl_vector_scale and memcpy.
    new_q->m_dofs[i] = current->m_dofs[i] + min(m_stepSize, currNorm) * gsl_vector_get(projected_gradient, i);
    //new_q->m_sumProjSteps[i] =
    //    min(m_stepSize, currNorm) * gsl_vector_get(projected_gradient, i) + current->m_sumProjSteps[i];
  }

  gsl_vector_free(projected_gradient);
  return new_q;
}
