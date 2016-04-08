//
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_GLOBALROTATEDOF_H
#define KGS_GLOBALROTATEDOF_H


#include "DOF.h"

class GlobalRotateDOF: public DOF {
 public:
  GlobalRotateDOF(const KinEdge* edge, int axis);

  Math3D::Vector3 getDerivative(Coordinate& coord) const;

  double getGlobalValue() const;

 protected:

  void updateEndVertexTransformation();

 private:
  int m_axis;
};


#endif //KGS_GLOBALROTATEDOF_H
