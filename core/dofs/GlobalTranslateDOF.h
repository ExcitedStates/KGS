//
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_GLOBALTRANSLATEDOF_H
#define KGS_GLOBALTRANSLATEDOF_H


#include "DOF.h"

class GlobalTranslateDOF: public DOF {
 public:
  GlobalTranslateDOF(const KinEdge* edge, int axis);

  Math3D::Vector3 getDerivative(Coordinate& coord) const override;

  double getGlobalValue() const override;

  double getRandomPerturbation() const override;

  double getMaxPerturbation() const override;

 protected:

  void updateEndVertexTransformation() override;

 private:
  int m_axis;
  static const double m_maxValue;
};


#endif //KGS_GLOBALTRANSLATEDOF_H
