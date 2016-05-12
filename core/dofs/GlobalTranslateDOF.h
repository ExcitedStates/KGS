//
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_GLOBALTRANSLATEDOF_H
#define KGS_GLOBALTRANSLATEDOF_H


#include "DOF.h"

class GlobalTranslateDOF: public DOF {
 public:
  GlobalTranslateDOF(const KinEdge* edge, int axis);

  Math3D::Vector3 getDerivative(Coordinate& coord) const;

  double getGlobalValue() const;

  double getMaxValue() const;
 protected:

  void updateEndVertexTransformation();

 private:
  int m_axis;
};


#endif //KGS_GLOBALTRANSLATEDOF_H
