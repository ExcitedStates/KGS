//
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_TORSIONDOF_H
#define KGS_TORSIONDOF_H


#include "DOF.h"

class TorsionDOF: public DOF {
 public:
  TorsionDOF(const KinEdge* edge): DOF(edge){}

  Math3D::Vector3 getDerivative(Coordinate& coord) const;

  double getGlobalValue() const;

  double getMaxValue() const;
 protected:

  void updateEndVertexTransformation();
};


#endif //KGS_TORSIONDOF_H
