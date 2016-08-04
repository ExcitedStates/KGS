#ifndef KGS_MKLQR_H
#define KGS_MKLQR_H

#include "math/QR.h"

class MKLQR: public QR {
public:
    MKLQR(gsl_matrix* M): QR(M) {}

    void updateFromMatrix() override;
};


#endif //KGS_MKLQR_H
