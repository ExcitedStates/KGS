//
// Created by Dominik Budday on 15.03.16.
//

#ifndef KGS_CLASHAVOIDINGMOVE_H
#define KGS_CLASHAVOIDINGMOVE_H

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "metrics/Metric.h"
#include "metrics/RMSD.h"
#include "moves/Move.h"
#include "core/Configuration.h"

class ClashAvoidingMove : public Move
{
 public:
  ClashAvoidingMove();

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:
  gsl_vector* projectOnClashNullspace(Configuration *conf,
                                      gsl_vector *gradient,
                                      std::set<std::pair<Atom *, Atom *> > &collisions);

  gsl_matrix* computeClashAvoidingJacobian( Configuration* conf,
                                            std::map<int,int>& dofMap,
                                            std::set< std::pair<Atom*,Atom*> >& collisions);

  /** Return a map that associates cycle-dofs and constrained dofs with a general dofs. */
  std::map<int,int> collectConstrainedDofMap(Configuration* conf, std::set< std::pair<Atom*,Atom*> >& allCollisions);
  const double m_maxRotation;
  const int m_trialSteps;
  const bool m_projectConstraints;
  const std::string m_collisionCheckAtomTypes;
};


#endif //KGS_CLASHAVOIDINGMOVE_H
