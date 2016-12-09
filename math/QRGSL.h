
#ifndef KGS_GSLQR_H
#define KGS_GSLQR_H

#include "math/QR.h"
#include <gsl/gsl_matrix.h>

class QRGSL: public QR {
 public:
  QRGSL(gsl_matrix* M): QR(M){}

 protected:

  void updateFromMatrix() override;
};


#endif //KGS_GSLQR_H
