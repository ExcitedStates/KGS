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
  gsl_matrix* computeClashAvoidingJacobian(Configuration* conf, std::set< std::pair<Atom*,Atom*> >& allCollisions, bool projectConstraints);
  const double m_maxRotation;
  const int m_trialSteps;
  const double m_stepSize;
  const std::string m_collisionCheck;
};


#endif //KGS_CLASHAVOIDINGMOVE_H
