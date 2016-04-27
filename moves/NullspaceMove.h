#ifndef NULLSPACEMOVE_H_
#define NULLSPACEMOVE_H_

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "metrics/Metric.h"
#include "metrics/RMSD.h"
#include "moves/Move.h"
#include "core/Configuration.h"


class NullspaceMove: public Move
{
 public:
  NullspaceMove(double maxRotation);

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:
  const double m_maxRotation;
  //const int m_decreaseSteps;
  //const int m_decreaseFactor;
  //const double m_stepSize;
};
#endif
