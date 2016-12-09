//
// Created by Rasmus Fonseca on 05/04/16.
//

#ifndef KGS_GLOBALROTATEDOF_H
#define KGS_GLOBALROTATEDOF_H


#include "DOF.h"

class GlobalRotateDOF: public DOF {
 public:
  GlobalRotateDOF(const KinEdge* edge, int axis);

  Math3D::Vector3 getDerivative(Coordinate& coord) const override;

  double getGlobalValue() const override;

  double getRandomPerturbation() const override;

  double getMaxPerturbation() const override;

 protected:

  void updateEndVertexTransformation() override;

 private:
  int m_axis;
  Atom* m_firstAtom;
  static const double m_maxValue;
};


#endif //KGS_GLOBALROTATEDOF_H
