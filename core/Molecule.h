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
#ifndef PROTEIN_H
#define PROTEIN_H

#include <map>
#include <string>
#include <list>
#include <set>
#include <core/graph/KinTree.h>

#include "Rigidbody.h"
#include "core/graph/KinGraph.h"
#include "core/Configuration.h"

class Chain;
class Grid;
class Bond;
class Hbond;

class Molecule {
 public:
  Molecule();
  ~Molecule();
  void setName(std::string& name);
  std::string getName() const;

  Atom* addAtom(const std::string& chain_name,
                const std::string& res_name,
                const int& res_id,
                const std::string& atomName,
                const int& atomId,
                const Coordinate& position );

  Chain* getChain (const std::string& chainName);
  Atom* getAtom (int atom_id);
  Atom* getAtom(const std::string& chainName, const int& resNum, const std::string& name);
  const std::vector<Atom*>& getAtoms() const;
  std::vector<Atom*>& getAtoms();
  const std::list<Bond*>& getCovBonds() const;
  std::list<Bond*>& getCovBonds();
  const std::list<Hbond*>& getHBonds() const;
  std::list<Hbond*>& getHBonds();

  int getMinResidueNumber();
  int getMaxResidueNumber();
  void printAllAtoms () const;
  int size() const;
  int totalDofNum () const;
  void printSummaryInfo() const;
  void updateAtom (int atom_id, Coordinate new_pos);
  bool inCollision (std::string collisionCheckAtoms = "all" ) const;
  std::set< std::pair<Atom*,Atom*> > getAllCollisions (std::string collisionCheckAtoms = "all" ) const;
  double minCollisionFactor (std::string collisionCheckAtoms = "all" ) const;
  void printAllCollisions () const;
  bool hasCycle() const;
//  void backupAtomIndex();
//  void restoreAtomIndex ();
  void alignReferencePositionsTo(Molecule * base);
  void translateReferencePositionsToRoot(Molecule * base);
  Grid* getGrid();
  void setCollisionFactor(double collisionFactor);

  void setConfiguration(Configuration *q);
  int countOriginalDofs () const;
  Coordinate centerOfMass () const;
  Coordinate centerOfGeometry () const;
  double checkCycleClosure(Configuration *q);//Todo: Move this to configuration, this is conf dependent, not topology

  void addCovBond (Bond * bond);
  void addHbond (Hbond * hb);
  void setToHbondIntersection (Molecule * p2);
  void buildSpanningTree();
  unsigned int findBestRigidBodyMatch(int rootRBId, Molecule * target = nullptr);
  void computeAtomJacobian (Atom* atom, gsl_matrix** jacobian);
//  gsl_vector* getEndEffectors();
  //void ProjectOnCycleNullSpace (gsl_vector *to_project, gsl_vector *after_project);


//  gsl_vector* vdwGradient ();
  std::pair<double,double> vdwEnergy (std::set< std::pair<Atom*,Atom*> >* allCollisions, std::string collisionCheck);
  double vdwEnergy (std::string collisionCheck);//compute vdw energy


  Configuration* resampleSugars(int startRes, int endRes, Configuration* cur, int aggression);
  Configuration* localRebuild(std::vector<int>& resetDOFs, std::vector<double>& resetValues, std::vector<int>& recloseDOFs, std::vector<int>& ignoreDOFs, Configuration* cur);

  std::vector<Chain*> chains;
  std::map<unsigned int,Rigidbody*> Rigidbody_map_by_id;
  std::set< std::pair<Atom*,Atom*> > m_initialCollisions; // collisions in the initial conformation stored in pairs of atoms, and use the smaller atom id as key.

  // Topology of rigid bodies
  KinTree *m_spanning_tree;
  Math3D::RigidTransform *m_Transformation; // cache: store the m_transformation of each rigid body group

  // Configuration
  Configuration *m_conf;
  Configuration *m_conf_backup;

  // Jacobian matrices containing all DOFs for updating atom positions
  gsl_matrix* AtomJacobian1;
  gsl_matrix* AtomJacobian2;
  gsl_matrix* AtomJacobian3;


  int* residueAnnotations;


 private:
  std::string name_;
  void _SetConfiguration(Configuration *q); // set the positions of atoms at configuration q (according to the spanning tree)
  void _SetConfiguration(Configuration *q, KinVertex* root, std::vector<KinVertex*>& subVerts);

  Chain* addChain (const std::string& chainName);

  void restoreAtomPos();

  void indexAtoms();
  Grid *m_grid;
//  Grid *m_backupGrid;
  std::list<Bond *> m_covBonds;
  std::list<Hbond *> m_hBonds;
  std::vector<Atom*> m_atoms;

  double m_collisionFactor;
};

#endif
