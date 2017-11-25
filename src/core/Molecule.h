/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

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
#include "DBond.h"

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
  std::vector<Chain*> m_chains; ///< Chains in the molecule. Each chain contains residues which in turn contains atoms

  Molecule();
  ~Molecule();

  void setName(const std::string& name);
  std::string getName() const;
  Chain* getChain (const std::string& chainName) const;
  const std::vector<Chain*>  getChains () const { return m_chains;};
  Atom* getAtom (int atom_id) const;
  Atom* getAtom(const std::string& chainName, const int& resNum, const std::string& name) const;
  const std::vector<Atom*>& getAtoms() const;
  std::vector<Atom*>& getAtoms();
  const std::list<Bond*>& getCovBonds() const;
  std::list<Bond*>& getCovBonds();
  const std::list<Hbond*>& getHBonds() const;
  const std::list<DBond*>& getDBonds() const;
  std::list<Hbond*>& getHBonds();
  std::list<DBond*>& getDBonds();

  int getMinResidueNumber();
  int getMaxResidueNumber();
  int size() const;
  int totalDofNum () const;
  bool inCollision (std::string collisionCheckAtoms = "all" );
  std::set< std::pair<Atom*,Atom*> > getAllCollisions (std::string collisionCheckAtoms = "all" );
  double minCollisionFactor (std::string collisionCheckAtoms = "all" );
  void printAllCollisions () ;
  double alignReferencePositionsTo(Molecule * base,Selection &sel);
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
  void addDBond (DBond * db);
  void setToHbondIntersection (Molecule * p2);
  std::vector<int> findBestRigidBodyMatch(std::vector<int> rootID, Molecule * target = nullptr);


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
  void buildRigidBodies (Selection& movingResidues, int collapseLevel = 1);
  void initializeTree(Selection& movingResidues,double collisionFactor = 1.0, const std::vector<int> &roots = {},Molecule* target = nullptr);

 private:
  std::string m_name;
  std::set< std::pair<Atom*,Atom*> > m_initialCollisions; ///< Colliding atom-pairs in the initial conformation
  Grid *m_grid;
  std::list<Bond *> m_covBonds;
  std::list<Hbond *> m_hBonds;
  std::list<DBond *> m_dBonds;
  std::vector<Atom*> m_atoms;
  std::map<unsigned int,Rigidbody*> m_rigidBodyMap; ///< Map for quickly looking up rigid bodies by id
  double m_collisionFactor;

  void _SetConfiguration(Configuration *q); // set the positions of atoms at configuration q (according to the spanning tree)
  void _SetConfiguration(Configuration *q, KinVertex* root, std::vector<KinVertex*>& subVerts);

  Chain* addChain (const std::string& chainName);
  Atom* addAtom(
      const bool& hetatm,
      const std::string& chain_name,
      const std::string& res_name,
      const int& res_id,
      const std::string& atomName,
      const int& atomId,
      const Coordinate& position
  );
  Bond* addCovBond(Atom* atom1, Atom* atom2);

  void sortHbonds();

  void restoreAtomPos();
  void indexAtoms();

  void buildSpanningTree(const std::vector<int>& rootIds);

  friend class IO; //The IO::readPdb function uses the addAtom function and sort hBonds function

};

#endif
