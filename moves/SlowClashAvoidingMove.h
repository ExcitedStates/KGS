//
// Created by Dominik Budday on 15.03.16.
//

#ifndef KGS_SLOWCLASHAVOIDINGMOVE_H
#define KGS_SLOWCLASHAVOIDINGMOVE_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "metrics/Metric.h"
#include "metrics/RMSD.h"
#include "moves/Move.h"
#include "core/Configuration.h"

class SlowClashAvoidingMove : public Move
{
 public:
  SlowClashAvoidingMove();

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:
  Configuration* projectOnClashNullspace(Configuration *conf,
                                      gsl_vector *gradient,
                                      std::set<std::pair<Atom *, Atom *> > &collisions);

  gsl_matrix* computeClashAvoidingJacobian( Configuration* conf,
                                            std::set< std::pair<Atom*,Atom*> >& allCollisions);

  /** Return a map that associates cycle-dofs and constrained dofs with a general dofs. */
  const double m_maxRotation;
  const int m_trialSteps;
  const bool m_projectConstraints;
  const std::string m_collisionCheckAtomTypes;
};


#endif //KGS_CLASHAVOIDINGMOVE_H
