//
// Created by Dominik Budday on 03.08.16.
//

#ifndef KGS_FIXEDLINK_H
#define KGS_FIXEDLINK_H

#include "DOF.h"

class FixedLink: public DOF {
 public:
  FixedLink(const KinEdge* edge): DOF(edge){}

  Math3D::Vector3 getDerivative(Coordinate& coord) const;

  double getGlobalValue() const;

  double getMaxValue() const;

  double getRandomPerturbation() const;

  double getMaxPerturbation() const;

 protected:

  void updateEndVertexTransformation();

 private:
  static const double m_maxValue;
};


#endif //KGS_FIXEDLINK_H
