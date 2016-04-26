#include "NullspaceMove.h"

#include <math/gsl_helpers.h>
#include "Logger.h"

using namespace std;

NullspaceMove::NullspaceMove():
    m_maxRotation(SamplingOptions::getOptions()->maxRotation),
    m_decreaseSteps(SamplingOptions::getOptions()->trialSteps),
    m_decreaseFactor(SamplingOptions::getOptions()->decreaseFactor),
    m_stepSize(SamplingOptions::getOptions()->stepSize)
{
  m_movesAccepted = 0;
  m_movesRejected = 0;
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

  //controlMaxRotation(projected_gradient);
  gsl_vector_scale_max_component(projected_gradient, SamplingOptions::getOptions()->maxRotation);

  //If resulting structure is in collision try scaling down the gradient
  for (int reduceStep = 0; reduceStep <= m_decreaseSteps; reduceStep++) {
    Configuration *new_q = new Configuration(current);

    for (int i = 0; i < new_q->m_numDOFs; ++i) {
      new_q->m_dofs[i] = current->m_dofs[i] + min(m_stepSize, currNorm) * gsl_vector_get(projected_gradient, i);
      new_q->m_sumProjSteps[i] =
          min(m_stepSize, currNorm) * gsl_vector_get(projected_gradient, i) + current->m_sumProjSteps[i];
    }

    if (new_q->updatedMolecule()->inCollision()) {
      log("dominik") << "Rejected!" << endl;
      m_movesRejected++;

      if (m_decreaseSteps > 0) {//reduce step size, for clash prevention use clash-avoiding move
        log("dominik") << "Scaling down with decrease factor" << endl;
        gsl_vector_scale(projected_gradient, m_decreaseFactor);
      }
      delete new_q;

    } else {//collision free
      m_movesAccepted++;

      gsl_vector_free(projected_gradient);

      return new_q;
    }
  }//end steps

  gsl_vector_free(projected_gradient);

  return nullptr;
}
