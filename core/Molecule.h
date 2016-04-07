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

#include "Rigidbody.h"
#include "RigidbodyGraph.h"
#include "core/Configuration.h"

class Chain;
class Grid;
class Bond;
class Hbond;
class Configuration;

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

  int getMinResidueNumber();
	int getMaxResidueNumber();
	void printAllAtoms () const;
	int size() const;
	int totalDofNum () const;
	void printSummaryInfo() const;
	void printBackboneAngleAndLength (std::string length_file="", std::string angle_file="") const;
	void updateAtom (int atom_id, Coordinate new_pos);
	void indexAtoms();
	bool inCollision (std::string collisionCheckAtoms = "all" ) const;
	std::set< std::pair<Atom*,Atom*> > getAllCollisions (std::string collisionCheckAtoms = "all" ) const;
	double minCollisionFactor (std::string collisionCheckAtoms = "all" ) const;
	void printAllCollisions () const;
	bool hasCycle() const;
	void backupAtomIndex(); 
	void restoreAtomIndex ();
	void alignReferencePositionsTo(Molecule * base);
	
	void SetConfiguration(Configuration *q);
	int countOriginalDofs () const;
	Coordinate centerOfMass () const;
	Coordinate centerOfGeometry () const;
	void checkCycleClosure(Configuration *q);
	
	void addCovBond (Bond * bond);
	void addHbond (Hbond * hb);
	void setToHbondIntersection (Molecule * p2);
	void buildRigidbodyTree(unsigned int bestRigidBody = 0 , bool flexibleSugar=true);
  unsigned int findBestRigidBodyMatch(int rootRBId, Molecule * target = NULL);
	RigidbodyGraphVertex* getRigidbodyGraphVertex (Atom* atom) const; // Return the vertex which is associated with the smallest DOF id edge among all the atom's rigidbodies.
	void computeAtomJacobian (Atom* atom, gsl_matrix** jacobian);
  gsl_vector* getEndEffectors();
	void ProjectOnCycleNullSpace (gsl_vector *to_project, gsl_vector *after_project);


	gsl_vector* vdwGradient ();
	std::pair<double,double> vdwEnergy (std::set< std::pair<Atom*,Atom*> >* allCollisions, std::string collisionCheck);

	Configuration* resampleSugars(int startRes, int endRes, Configuration* cur, int aggression);
	Configuration* localRebuild(std::vector<int>& resetDOFs, std::vector<double>& resetValues, std::vector<int>& recloseDOFs, std::vector<int>& ignoreDOFs, Configuration* cur);

	std::vector<Atom*> atoms;
	std::vector<Chain*> chains;
	std::map<unsigned int,Rigidbody*> Rigidbody_map_by_id;
	Grid *Atom_pos_index;
	Grid *backup_Atom_pos_index;
	std::list<Bond *> Cov_bonds;
	std::list<Hbond *> H_bonds; // To do: or it is better to use list<>?
	std::set< std::pair<Atom*,Atom*> > Initial_collisions; // collisions in the initial conformation stored in pairs of atoms, and use the smaller atom id as key.

	// Topology of rigid bodies
	RigidbodyTree *m_spanning_tree;
	RigidTransform *m_Transformation; // cache: store the transformation of each rigid body group

	// Configuration
	Configuration *m_conf;
	Configuration *m_conf_backup;

	// Jabocian matrices containing all DOFs for updating atom positions
	gsl_matrix* AtomJacobian1;
	gsl_matrix* AtomJacobian2;
	gsl_matrix* AtomJacobian3;


	int* residueAnnotations;


private:
	std::string name_;
	void _SetConfiguration(Configuration *q); // set the positions of atoms at configuration q (according to the spanning tree)
	void _SetConfiguration(Configuration *q, RigidbodyGraphVertex* root, std::vector<RigidbodyGraphVertex*>& subVerts);

  Chain* addChain (const std::string& chainName);

  void RestoreAtomPos();
};

#endif
