#include "Move.h"
#include <math/gsl_helpers.h>
#include "Logger.h"


Move::~Move(){}

Configuration* Move::move(Configuration* current, gsl_vector* gradient)
{
  return performMove(current, gradient);
}

