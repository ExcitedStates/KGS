#ifndef COMPOSITE_MOVE_
#define COMPOSITE_MOVE_

#include <vector>
#include <utility>

#include "Move.h"

#include "core/Configuration.h"

class CompositeMove: public Move
{
 public:
  CompositeMove();
  ~CompositeMove();

  void addMove(Move* m, double weight);

 protected:
  Configuration* performMove(Configuration* current, gsl_vector* gradient);

 private:

  std::vector< std::pair<Move*, double> > child_moves;

  double total_weight;
};

#endif
