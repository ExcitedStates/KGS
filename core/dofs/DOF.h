
#ifndef KGS_DOF_H
#define KGS_DOF_H


#include <core/graph/KinGraph.h>

/**
 * A degree-of-freedom affecting a branch of the RigidbodyTree.
 * All DOFs must be able to return an absolute value (for instance the N-CA-C-N dihedral angle) but in practice
 * relative angles (related to some conformational origo) are stored and updated when working with DOFs.
 *
 * Each DOF is responsible for returning derivatives of coordinates 'downstream' in the kinematic tree. Additionally,
 * they are responsible for forward-propagation of transformations.
 */
class DOF {
 public:
  DOF(Edge* edge);

  /**
   * Return the partial derivative of the specified position assuming that `coord` is affected
   * by this DOF. No checks are performed to verify this.
   */
  virtual Math3D::Vector3 getDerivative(Coordinate& coord) const = 0;

  /**
   * Set the relative value of this DOF.
   * Calling this function "schedules" the change but nothing substantial is actually performed before
   * `forwardPropagate` is called.
   */
  void setValue(double delta);

  /**
   * Get the relative value of this DOF
   */
  double getValue() const;

  /**
   * Get the global value of this DOF.
   */
  virtual double getGlobalValue() const = 0;

 protected:
//  /**
//   * Update the transformation `tr`, apply it to `m_edge->EndVertex`, and propagate to edges leaving
//   * `m_edge->EndVertex`.
//   */
//  virtual void forwardPropagate(Math3D::RigidTransform& tr)=0;

  /**
   * Update the transformation matrix in `m_edge->EndVertex` based on the one in `m_edge->StartVertex` and
   * the value for this DOF.
   */
  virtual void updateEndVertexTransformation() = 0;

  Edge const * m_edge;

  double m_value;

};


#endif //KGS_DOF_H
