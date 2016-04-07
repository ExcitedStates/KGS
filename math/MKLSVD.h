#ifndef KGS_MKLSVD_H
#define KGS_MKLSVD_H

#include "math/SVD.h"

class MKLSVD: public SVD {
public:
    MKLSVD(gsl_matrix* M): SVD(M) {}

    void UpdateFromMatrix() override;
};


#endif //KGS_MKLSVD_H
