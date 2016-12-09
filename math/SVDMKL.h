#ifndef KGS_SVDMKL_H
#define KGS_SVDMKL_H

#include "math/SVD.h"

class SVDMKL: public SVD {
public:
    SVDMKL(gsl_matrix* M): SVD(M) {}

    void UpdateFromMatrix() override;
};


#endif //KGS_SVDMKL_H
