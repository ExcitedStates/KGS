/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

        Permission is hereby granted, free of charge, to any person obtaining a copy of
        this software and associated documentation files (the "Software"), to deal in
        the Software without restriction, including without limitation the rights to
        use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
        of the Software, and to permit persons to whom the Software is furnished to do
        so, subject to the following conditions:

        This entire text, including the above copyright notice and this permission notice
        shall be included in all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
        OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
        FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
        IN THE SOFTWARE.


*/
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <list>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <vector>
#include <math/Nullspace.h>

#include "core/graph/KinGraph.h"
#include "core/graph/SugarVertex.h"

class Molecule;
class Configuration;

typedef std::list<Configuration*> ConfigurationList;
typedef std::vector<Configuration*> ConfigurationArray;


/**
 * A configuration holds all DOF-values necessary to update a molecules atom-positions.
 * DOF-values are relative to whatever is stored in the Atom::reference_positions, so calling
 * Configuration* conf = new Configuration(m_protein);
 * will create a configuration where all DOF-values are 0 and which represents whatever is in
 * the reference positions of m_protein. Modifying for example
 * conf->m_dofs[2] += 0.1
 * will add 0.1 units (often radians) to the third DOF. To see the resulting structure call either
 * m_protein->SetConfiguration(conf);
 * or
 * conf->updatedProtein();
 * and then access Atom::position (not Atom::m_referencePosition).
 *
 * As new configurations are generated from older ones, the field m_children and m_parent store the
 * connectivity information.
 */
class Configuration
{
 public:
  int m_numDOFs;                    ///< Number of DOFs
  double *m_dofs;                    ///< DOF-values (relative to Atom::m_referencePosition)
  double *m_sumProjSteps;          //TODO: What is this?

  /** Construct a configuration with all DOF-values set to 0 and no m_parent. */
  Configuration(Molecule * protein);

  /** Construct a configuration with all DOF-values set to 0 and the specified m_parent. */
  Configuration(Configuration* parent);

  ~Configuration();

  void updateGlobalTorsions();           ///< Update the global DOF-values (m_dofs_global field)
  double getGlobalTorsions(int i) const; ///< Get a global DOF-value
  double* getGlobalTorsions() const;     ///< Get global DOF-value array

  /** Set the specified dof to the global torsion value. Convenient function for
   * calculating difference between reference torsion and val. */
  void setGlobalTorsion(int i, double val);

  Configuration* clone() const;          ///< Copy this configuration
  void SetAll(double v);                 // TODO: Remove (?)
  void Normalize();                      // TODO: Remove (?)
  double Length();                       // TODO: Remove (?)

  void Print();                          // TODO: Remove or rename to printDOFs
  void computeCycleJacobianAndNullSpace();
  void computeJacobians();               ///< Compute non-redundant cycle jacobian and hbond-jacobian
  void identifyBiggerRigidBodies();      ///< Identify clusters
  void readBiggerSet();                  ///< read the set of clusters, related to identifying clusters

  //Clash-avoiding Jacobian and nullspace
  void ComputeClashAvoidingJacobianAndNullSpace (std::map< std::pair<Atom*,Atom*>,int >  allCollisions,bool firstTime,bool projectConstraints);
  void computeClashAvoidingJacobian (std::map< std::pair<Atom*,Atom*>,int > allCollisions,bool projectConstraints);


  static bool compareSize(std::pair<int, unsigned int> firstEntry, std::pair<int, unsigned int> secondEntry);//TODO: What is this?


  // When the samples are generated as an expanding tree, the m_children and m_parent store the connectivity information of these nodes
  const int m_treeDepth;             ///< Depth in the exploration tree
  int m_id;                         ///< ID of configuration
  double m_vdwEnergy;               ///< van der Waals energy of configuration
  double m_deltaH;                  ///< change in enthalpy due to clash constraints
  double m_distanceToTarget;        ///< Distance to target configuration
  double m_distanceToParent;        ///< Distance to m_parent configuration
  double m_distanceToIni;           ///< Distance to initial configuration
  double m_paretoFrontDistance;
  double m_maxConstraintViolation;  //Maximum detected distance violation of h-bond constraints, Todo: maybe not necessary to keep
  double m_minCollisionFactor;      //minimum necessary clash-factor for configuration to be clash free, Todo: maybe not necessary to keep
  double m_usedClashPrevention;

  std::map<unsigned int, Rigidbody*> m_biggerRBMap;  // TODO: What is the int? Cluster-idx?
  std::vector< std::pair<int, unsigned int> > m_sortedRBs; // TODO: Sorted on what? What does the integers refer to?

  int m_numClusters;            ///< Number of rigid clusters (super rigid bodies)
  int m_maxIndex;               ///< Index of largest cluster
  int m_maxSize;                ///< Size of largest cluster
  int m_clashFreeDofs;          ///< Number of clash-free dofs (for post-processing)

  Molecule * updatedProtein();    ///< Update the atom-positions to reflect this configuration and return the m_protein
  Molecule * getProtein() const;  ///< Return the associated m_protein
  void updateProtein();         ///< Update the atom-positions to reflect this configuration

  /** Return the cycle jacobian. Assumes that computeJacobians has been called on this configuration last */
  gsl_matrix* getCycleJacobian() const;

  Nullspace* getNullspace();    ///< Compute the nullspace (if it wasn't already) and return it

  Configuration* getParent();   ///< Access configuration that spawned this one
  ConfigurationList& getChildren(); ///< Access child configurations

  static Nullspace* ClashAvoidingNullSpace; //TODO: Make private (or even better put in ClashAvoidingMove).
 protected:

  double *m_dofs_global;                 ///< DOF-values in a global system (not relative to Atom::reference_position)
  Molecule * const m_protein;            ///< The m_protein related to the configuration
  Configuration * const m_parent;        ///< The m_parent-configuration this configuration was generated from
  ConfigurationList m_children;          ///< List of child-configurations

  // Jacobian matrix of all the cycles of rigid bodies
  static gsl_matrix* CycleJacobian; // column dimension is the number of DOFs; row dimension is 5 times the number of cycles because 2 atoms on each cycle-closing edge
  static gsl_matrix* HBondJacobian; // column dimension is the number of DOFS; row dimension is the number of cycles
  static gsl_matrix* ClashAvoidingJacobian;

  static SVD* JacobianSVD;
  Nullspace* nullspace;                  ///< Nullspace of this configuration
};




#endif

