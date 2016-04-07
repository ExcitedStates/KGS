
#ifndef KGS_GSLSVD_H
#define KGS_GSLSVD_H

#include "math/SVD.h"
#include <gsl/gsl_matrix.h>

class GSLSVD: public SVD {
 public:
  GSLSVD(gsl_matrix* M): SVD(M){}

 protected:

  void UpdateFromMatrix() override;
};


#endif //KGS_GSLSVD_H
