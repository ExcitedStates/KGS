#ifndef KGS_CUDASVD_H
#define KGS_CUDASVD_H

#include "math/SVD.h"

class CudaSVD: public SVD {
 public:
  CudaSVD(gsl_matrix* M): SVD(M) {}

  void UpdateFromMatrix() override;
 private:
};


#endif //KGS_CUDASVD_H
