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

#ifndef GRID_H
#define GRID_H

#include <map>
#include <vector>

#include "Molecule.h"

class Molecule;


#define GRID_CELL_SIZE 4.0
//extern double COLLISION_FACTOR;
extern double RADIUS_RATIO;

struct Coordinate_cmp {
	bool operator() (const Coordinate& c1, const Coordinate& c2) const {
		//New version: 5 comparisons worst-case, 1.5 average case @RFonseca
		if (c1.x<c2.x)	return true;
		if (c1.x>c2.x)	return false;
		if (c1.y<c2.y)	return true;
		if (c1.y>c2.y)	return false;
		return c1.z<c2.z;

		//Original version: 7 comparisons worst-case, 3.5 average-case @RFonseca
		if (c1.x<c2.x)	return true;
		if (c1.x>c2.x)	return false;
		else if (c1.x==c2.x && c1.y<c2.y)
			return true;
		else if (c1.x==c2.x && c1.y==c2.y && c1.z<c2.z)
			return true;
		return false;
	}
};

class Grid {
 public:
  Grid (Molecule * protein, double collisionFactor=1.0);
  Grid ();
  ~Grid ();
  void print() const;
  std::vector<Atom*> getNeighboringAtoms (Atom* atom, bool neighborWithLargerId=true, bool noCovBondNeighbor=true, bool noHbondNeighbor=true, double radius=GRID_CELL_SIZE) const;
  std::vector<Atom*> getNeighboringAtomsVDW (Atom* atom, bool neighborWithLargerId=true, bool noCovBondNeighbor=true, bool noSecondCovBondNeighbor=true, bool noHbondNeighbor=true, double radius=GRID_CELL_SIZE) const;

  bool inCollision (Atom* atom, std::set< std::pair<Atom*,Atom*> > const &initial_collision_list, std::string collisionCheckAtoms= "all", bool onlyCheckLargerIds=true) const;
  double minFactorWithoutCollision (Atom* atom, std::set< std::pair<Atom*,Atom*> > const &initial_collision_list, std::string collisionCheckAtoms= "all", bool onlyCheckLargerIds=true) const;
  std::vector<Atom*> getAllCollisions (Atom* atom, std::set< std::pair<Atom*,Atom*> > const &initial_collision_list, std::string collisionCheckAtoms= "all",bool onlyCheckLargerIds=true) const;

  bool removeAtom (Atom* atom);
  void addAtom (Atom* atom);

	void setCollisionFactor(double collisionFactor);

 private:
	Coordinate makeKey (Coordinate& pos) const;
  static const double Cell_size;
  double m_collisionFactor;

  double Max_x;
  double Min_x;
  double Max_y;
  double Min_y;
  double Max_z;
  double Min_z;

  std::map<Coordinate,std::vector<Atom*>*,Coordinate_cmp> Atom_map;
};

#endif
