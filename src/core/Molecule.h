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
#include <Selection.h>

#include "Rigidbody.h"
#include "core/graph/KinGraph.h"
#include "core/Configuration.h"

class Chain;
class Grid;
class Bond;
class Hbond;
class Selection;

class Molecule {
 public:

  KinTree *m_spanningTree; ///< Topology of rigid bodies
  Configuration *m_conf;   ///< Currently set configuration
  int* residueAnnotations;
  std::vector<Chain*> chains; ///< Chains in the molecule. Each chain contains residues which in turn contains atoms

  Molecule();
  ~Molecule();

  void setName(const std::string& name);
  std::string getName() const;
  Chain* getChain (const std::string& chainName) const;
  Atom* getAtom (int atom_id) const;
  Atom* getAtom(const std::string& chainName, const int& resNum, const std::string& name) const;
  const std::vector<Atom*>& getAtoms() const;
  std::vector<Atom*>& getAtoms();
  const std::list<Bond*>& getCovBonds() const;
  std::list<Bond*>& getCovBonds();
  const std::list<Hbond*>& getHBonds() const;
  std::list<Hbond*>& getHBonds();

  int getMinResidueNumber();
  int getMaxResidueNumber();
  int size() const;
  int totalDofNum () const;
  bool inCollision (std::string collisionCheckAtoms = "all" ) const;
  std::set< std::pair<Atom*,Atom*> > getAllCollisions (std::string collisionCheckAtoms = "all" ) const;
  double minCollisionFactor (std::string collisionCheckAtoms = "all" ) const;
  void printAllCollisions () const;
  void alignReferencePositionsTo(Molecule * base);
  void translateReferencePositionsToRoot(Molecule * base);
  Grid* getGrid();
  void setCollisionFactor(double collisionFactor);

  void forceUpdateConfiguration(Configuration *q);
  void setConfiguration(Configuration *q);
  int countOriginalDofs () const;
  Coordinate centerOfMass () const;
  Coordinate centerOfGeometry () const;
  double checkCycleClosure(Configuration *q);//Todo: Move this to configuration, this is conf dependent, not topology
  void computeCycleViolation(Configuration *q, gsl_vector *currentViolation);

  void addCovBond (Bond * bond);
  void addHbond (Hbond * hb);
  void setToHbondIntersection (Molecule * p2);
  unsigned int findBestRigidBodyMatch(int rootRBId, Molecule * target = nullptr);


  std::pair<double,double> vdwEnergy (std::set< std::pair<Atom*,Atom*> >* allCollisions, std::string collisionCheck);
  double vdwEnergy (std::string collisionCheck);//compute vdw energy

  std::set< std::pair<Atom*,Atom*> >& getInitialCollisions(); ///< Colliding atom-pairs in the initial conformation

  Configuration* resampleSugars(int startRes, int endRes, Configuration* cur, int aggression);
  Configuration* localRebuild(std::vector<int>& resetDOFs, std::vector<double>& resetValues, std::vector<int>& recloseDOFs, std::vector<int>& ignoreDOFs, Configuration* cur);

  /** Clone atom positions and bonds, but don't create rigid bodies and spanning tree. */
  Molecule* deepClone() const;

  Molecule* collapseRigidBonds(int collapseLevel = 1);

  void writeRigidbodyIDToBFactor();

  const std::vector<Rigidbody*> getRigidbodies() const;
  void buildRigidBodies (Selection& movingResidues);

 private:
  std::string m_name;
  std::set< std::pair<Atom*,Atom*> > m_initialCollisions; ///< Colliding atom-pairs in the initial conformation
  Grid *m_grid;
  std::list<Bond *> m_covBonds;
  std::list<Hbond *> m_hBonds;
  std::vector<Atom*> m_atoms;
  std::map<unsigned int,Rigidbody*> m_rigidBodyMap; ///< Map for quickly looking up rigid bodies by id
  double m_collisionFactor;

  void _SetConfiguration(Configuration *q); // set the positions of atoms at configuration q (according to the spanning tree)
  void _SetConfiguration(Configuration *q, KinVertex* root, std::vector<KinVertex*>& subVerts);

  Chain* addChain (const std::string& chainName);
  Atom* addAtom(const std::string& chain_name,
                const std::string& res_name,
                const int& res_id,
                const std::string& atomName,
                const int& atomId,
                const Coordinate& position );
  Bond* addCovBond(Atom* atom1, Atom* atom2);

  void restoreAtomPos();
  void indexAtoms();

  void buildSpanningTree(const std::vector<int>& rootIds);

  friend class IO; //The IO::readPdb function uses the addAtom function

};

#endif
