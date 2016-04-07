
#ifndef REBUILDMOVE_H_
#define REBUILDMOVE_H_

#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "moves/Move.h"
#include "core/Configuration.h"
#include "SamplingOptions.h"

class Resampler{
 public:
//  virtual
};

class RebuildMove: public Move
{
 public:
  RebuildMove();

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:
  const int m_fragmentLength;
  const int m_rebuildAggression;

};
#endif
