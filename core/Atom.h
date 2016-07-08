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
#ifndef ATOM_H
#define ATOM_H

//#pragma warning (disable : 4244)
//#pragma warning (disable : 4267)

#include <string>
#include <vector>
#include "Residue.h"
#include "Coordinate.h"

// Vdw radii are from source http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html#allatom, except for SE and UNKNOWN
// To get CHARMM parameters for each atom; for now assume these VDW_RADII = Rmin/2
#define VDW_RADIUS_C 1.700
#define VDW_RADIUS_H 1.000
#define VDW_RADIUS_N 1.625
#define VDW_RADIUS_O 1.480
#define VDW_RADIUS_S 1.782
#define VDW_RADIUS_SE 1.880
#define VDW_RADIUS_UNKNOWN 2.000
// Epsilon parameters for each atom are taken from CHARMM FF by averaging all atom types per element
#define VDW_EPSILON_C -0.0661154
#define VDW_EPSILON_H -0.0361867 
#define VDW_EPSILON_N -0.2  
#define VDW_EPSILON_O -0.137663  
#define VDW_EPSILON_S -0.433333  
#define VDW_EPSILON_SE -0.433333 
#define VDW_EPSILON_UNKNOWN -0.2


// Mass are from source http://www.okstate.edu/jgelder/PTbyAtMass.GIF, the unknown mass is the max mass.
#define MASS_C 12.01
#define MASS_H 1.01
#define MASS_N 14.01
#define MASS_O 16.00
#define MASS_S 32.06
#define MASS_SE 78.96
#define MASS_UNKNOWN 267.00

typedef enum {
	atomC   = 0,
	atomH   = 1,
	atomN   = 2,
	atomO   = 3,
	atomS   = 4,
	atomSE  = 5,
	atomOther = 6,
	atomTypeAll = 7
} AtomType;

class Bond;
class Hbond;
class Rigidbody;


class Atom {
 public:

	Atom(const std::string& name, const int& id, const Coordinate& pos, Residue* residue);

	const std::string& getName() const;
	int getId() const;
	double getMass() const;    ///< Return the atomic mass (depends on element)
  double getRadius() const;  ///< Return the van der Waals radius (depends on element)
  double getEpsilon() const; ///< Return the van der Waals radius (depends on element)
  Residue* getResidue() const;

	void printSummaryInfo() const;


	void addCovBond (Bond * cbond);
  Atom* getIthCovNeighbor(int i) const;
  bool isCovNeighbor(Atom* other) const;
  bool isSecondCovNeighbor(Atom* other) const;
  Atom* getFirstCovNeighbor() const;

	void addHbond(Hbond * hbond);
	void removeHbond(Hbond * hbond);
  bool isHbondNeighbor(Atom* other) const;

	Atom* getBondNeighbor(Bond * bond) const;

  //TODO: Let the Coordinate class take care of this.
	double distanceTo(Atom* other) const; // Euclidean distance between atom self and other
  //TODO: Let the Coordinate class take care of this. This impl. is also slow (uses pow)
	bool isWithinDistanceFrom(Atom* center, double dist) const; // is within a sphere centered at center with radius dist

	void setRigidbody (Rigidbody* rb);
  Rigidbody* getRigidbody() const;
  void setBiggerRigidbody(Rigidbody* rb);
  Rigidbody* getBiggerRigidbody() const;
  bool inSameRigidbody(Atom* another) const;

	std::string getType() const; //TODO: Return type should be char or AtomType
  void setAsMainchainAtom();
  void setAsSidechainAtom();
  bool isSidechainAtom() const;
  bool isBackboneAtom() const;
  bool isHeavyAtom() const;
  bool isCollisionCheckAtom (std::string collisionCheckAtoms="all") const;

	bool compareName(std::string atomName) const;
	bool compareType(std::string atomType) const;
  static bool compare(Atom* atom1, Atom* atom2);

	Coordinate m_position;
	Coordinate m_referencePosition;
	bool On_sidechain;
	std::vector<Bond *> Cov_bond_list;
	std::vector<Hbond *> Hbond_list;
	std::vector<Atom*> Cov_neighbor_list;
	std::vector<Atom*> Hbond_neighbor_list;
	std::vector<Atom*> Second_cov_neighbor_list; // 2nd immediate covalent bond neighbors

  AtomType Element;

 private:
  const int Id;
  const std::string Name;
	Residue* m_parentResidue;
  Rigidbody* m_rigidbody;
  Rigidbody* m_biggerRigidbody;
};

std::ostream& operator<<(std::ostream& os, const Atom& a);
std::ostream& operator<<(std::ostream& os, const Atom* a);

#endif
