#ifndef KGS_CULASVD_H
#define KGS_CULASVD_H

#include "math/SVD.h"

#ifdef __GPU_CULA
#include <cula.hpp>
#endif

class CULASVD: public SVD {
 public:
  CULASVD(gsl_matrix* M): SVD(M) {}

  void UpdateFromMatrix() override;
 private:
  void printCulaError(culaStatus status);
  int meetsMinimumCulaRequirements();
};


#endif //KGS_CULASVD_H
