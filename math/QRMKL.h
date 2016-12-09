#ifndef KGS_QRMKL_H
#define KGS_QRMKL_H

#include "math/QR.h"

class QRMKL: public QR {
public:
    QRMKL(gsl_matrix* M): QR(M) {}

    void updateFromMatrix() override;
};


#endif //KGS_MKLQR_H
