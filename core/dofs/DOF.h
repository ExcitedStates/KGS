
#ifndef KGS_DOF_H
#define KGS_DOF_H


#include <core/graph/KinGraph.h>
#include <core/graph/KinEdge.h>

class KinEdge;

/**
 * A degree-of-freedom affecting a branch of the KinTree.
 * All DOFs must be able to return an absolute value (for instance the N-CA-C-N dihedral angle) but in practice
 * relative angles (related to some conformational origo) are stored and updated when working with DOFs.
 *
 * Each DOF is responsible for returning derivatives of coordinates 'downstream' in the kinematic tree. Additionally,
 * they are responsible for forward-propagation of transformations.
 */
class DOF {
 public:
  DOF(const KinEdge* edge);

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

  virtual double getMaxValue() const = 0;

  /**
   * Get the global value of this DOF.
   */
  virtual double getGlobalValue() const = 0;

  unsigned int getIndex() const;

  unsigned int getCycleIndex() const;

 protected:
  friend class KinEdge;
  friend class KinTree;

  /**
   * Update the m_transformation matrix in `m_edge->EndVertex` based on the one in `m_edge->StartVertex` and
   * the value for this DOF.
   */
  virtual void updateEndVertexTransformation() = 0;

  /** Set the DOF index. Only called from KinTree. */
  void setIndex(unsigned int idx);

  /** Set the cycle-DOF index. Only called from KinTree. */
  void setCycleIndex(unsigned int idx);

//  /**
//   * Update the m_transformation `tr`, apply it to `m_edge->EndVertex`, and propagate to m_edges leaving
//   * `m_edge->EndVertex`.
//   */
//  virtual void forwardPropagate(Math3D::RigidTransform& tr)=0;


  KinEdge const * m_edge;

  double m_value;

  int m_index;
  int m_cycleIndex;

};


#endif //KGS_DOF_H
