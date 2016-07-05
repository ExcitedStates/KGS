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
#include <fstream>
#include <iostream>
#include <math.h>

#include "Grid.h"

using namespace std;

const double Grid::Cell_size = GRID_CELL_SIZE;

Grid::Grid (Molecule * protein, double collisionFactor):
    m_collisionFactor(collisionFactor)
{
	Max_x = -1000;
	Max_y = -1000;
	Max_z = -1000;
	Min_x = 1000;
	Min_y = 1000;
	Min_z = 1000;
	for (Atom* const& atom: protein->getAtoms()) {
		if (atom->m_Position.x<Min_x) Min_x = atom->m_Position.x;
		if (atom->m_Position.y<Min_y) Min_y = atom->m_Position.y;
		if (atom->m_Position.z<Min_z) Min_z = atom->m_Position.z;
		if (atom->m_Position.x>Max_x) Max_x = atom->m_Position.x;
		if (atom->m_Position.y>Max_y) Max_y = atom->m_Position.y;
		if (atom->m_Position.z>Max_z) Max_z = atom->m_Position.z;
	}
	//cout << "InfoGridConstructor)\t (Min_x,Max_x) (Min_y,Max_y) (Min_z,Max_z) = (" << Min_x << "," << Max_x << ") ("
	//									          << Min_y << "," << Max_y << ") ("
	//										  << Min_z << "," << Max_z << ") (" << endl;

	for (Atom* const& atom: protein->getAtoms()) {
		addAtom( atom );
	}
}

Grid::~Grid () {
	for (auto const& coord_vec_pair: Atom_map) {
		delete coord_vec_pair.second;
	}
}

Grid::Grid(): m_collisionFactor(1.0) {
}

Grid* Grid::deepClone () const {
	Grid* grid_copy = new Grid();

	grid_copy->Max_x = Max_x;
	grid_copy->Min_x = Min_x;
	grid_copy->Max_y = Max_y;
	grid_copy->Min_y = Min_y;
	grid_copy->Max_z = Max_z;
	grid_copy->Min_z = Min_z;

	for (map<Coordinate,vector<Atom*>*,Coordinate_cmp>::const_iterator mit=Atom_map.begin(); mit!=Atom_map.end(); ++mit) {
		vector<Atom*> *list = new vector<Atom*>();
		for (vector<Atom*>::const_iterator vit=mit->second->begin(); vit!=mit->second->end(); ++vit)
			list->push_back(*vit);
		grid_copy->Atom_map[mit->first] = list;
	}

	return grid_copy;
}

void Grid::print() const {
	for (map<Coordinate,vector<Atom*>*,Coordinate_cmp>::const_iterator mit=Atom_map.begin(); mit!=Atom_map.end(); ++mit) {
		cout << "Key (" << mit->first.tostring() << ") => ";
		for (vector<Atom*>::const_iterator vit=mit->second->begin(); vit!=mit->second->end(); ++vit) 
			cout << (*vit)->getId() << " ";
		cout << endl;
	}
}

Coordinate Grid::makeKey (Coordinate& pos) const {
	int x = floor((pos.x-Min_x)/Cell_size);
	int y = floor((pos.y-Min_y)/Cell_size);
	int z = floor((pos.z-Min_z)/Cell_size);
	Coordinate key(x,y,z);
	return key;
}
//---------------------------------------------------------
vector<Atom*> Grid::getNeighboringAtoms (Atom* atom, bool neighborWithLargerId, bool noCovBondNeighbor, bool noHbondNeighbor, double radius) const {
// TODO: can optimize this function by not visiting some corner cells which is further away than radius*sqrt(2)
	vector<Atom*> neighbors;
	Coordinate pos = atom->m_Position;
	Coordinate key = makeKey(pos);
  double radSq = radius*radius;
	int neighbor_cell_num = int(ceil(radius/Cell_size));
	for (int i=-neighbor_cell_num; i<=neighbor_cell_num; ++i) {
		for (int j=-neighbor_cell_num; j<=neighbor_cell_num; ++j) {
			for (int k=-neighbor_cell_num; k<=neighbor_cell_num; ++k) {
				Coordinate cur_key(key.x+i,key.y+j,key.z+k);
				map<Coordinate,vector<Atom*>*,Coordinate_cmp>::const_iterator map_itr = Atom_map.find(cur_key);
				if ( map_itr == Atom_map.end() ) continue; // if the key doesn't exist
				for (vector<Atom*>::const_iterator atom_itr=map_itr->second->begin(); atom_itr!=map_itr->second->end(); ++atom_itr) {
					if ( (*atom_itr) != atom ) {
						if (neighborWithLargerId && (*atom_itr)->getId()<atom->getId())
							continue;
						if (noCovBondNeighbor && atom->isCovNeighbor(*atom_itr))
							continue;
						if (noHbondNeighbor && atom->isHbondNeighbor(*atom_itr))
							continue;
//            double dist = atom->distanceTo(*atom_itr);
//            if (dist>radius) continue;
						if (atom->m_Position.distanceSquared( (*atom_itr)->m_Position )>radSq) continue;
						neighbors.push_back(*atom_itr);
					}
				}
			}
		}
	}
	return neighbors;
}
//---------------------------------------------------------
vector<Atom*> Grid::getNeighboringAtomsVDW (Atom* atom, bool neighborWithLargerId, bool noCovBondNeighbor, bool noSecondCovBondNeighbor, bool noHbondNeighbor, double radius) const {
// TODO: can optimize this function by not visiting some corner cells which is further away than radius*sqrt(2)
	vector<Atom*> neighbors;
	Coordinate pos = atom->m_Position;
	Coordinate key = makeKey(pos);
	int neighbor_cell_num = int(ceil(radius/Cell_size)); // radius = Neighbor list cutoff; defaults to Cell_Size
	for (int i=-neighbor_cell_num; i<=neighbor_cell_num; ++i)
		for (int j=-neighbor_cell_num; j<=neighbor_cell_num; ++j)
			for (int k=-neighbor_cell_num; k<=neighbor_cell_num; ++k) {
				Coordinate cur_key(key.x+i,key.y+j,key.z+k);
				map<Coordinate,vector<Atom*>*,Coordinate_cmp>::const_iterator map_itr = Atom_map.find(cur_key);
				if ( map_itr == Atom_map.end() ) continue; // if the key doesn't exist
				for (vector<Atom*>::const_iterator atom_itr=map_itr->second->begin(); atom_itr!=map_itr->second->end(); ++atom_itr)
//					if ( (*atom_itr) != atom && (*atom_itr)->isWithinDistanceFrom( atom, radius ) ) {
          if ( (*atom_itr) != atom && (*atom_itr)->m_Position.distanceSquared(atom->m_Position)<=radius*radius ) {
					//	if (neighborWithLargerId && (*atom_itr)->Id<atom->Id)
					//		continue;
						if (noCovBondNeighbor && atom->isCovNeighbor(*atom_itr))
							continue;
						if (noSecondCovBondNeighbor && atom->isSecondCovNeighbor(*atom_itr))
							continue;
					//	if (noHbondNeighbor && atom->isHbondNeighbor(*atom_itr))
					//		continue;
					//	double dist = atom->distanceTo(*atom_itr);
					//	if (dist>radius) continue;
						neighbors.push_back(*atom_itr);
					}
			}
	return neighbors;
}
//---------------------------------------------------------
/*
 * Given an atom and an initial list of collisions, determine if the atom is colliding with another atom
 * by ignoring atom pairs in the initial list.
 */
bool Grid::inCollision (Atom* atom, set< pair<Atom*,Atom*> > const &initial_collision_list, std::string collisionCheckAtoms, bool onlyCheckLargerIds) const {

	vector<Atom*> neighbors = getNeighboringAtoms(atom,onlyCheckLargerIds);
	//cout << "\t ------------------- Checking collisions for atom ... " << atom->getResidue()->getId() << " " <<  atom->getName() << " " <<  atom->Id;
	//cout << "   with " << neighbors.size() << " neighbors." << endl;
	for (vector<Atom*>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it) {
		if ( atom->isCovNeighbor(*it) || atom->isSecondCovNeighbor(*it) || atom->isHbondNeighbor(*it) || !(atom->isCollisionCheckAtom(collisionCheckAtoms))) {
			continue;
		}
		double dist = atom->m_Position.distanceTo((*it)->m_Position);
		double threshold =  m_collisionFactor * ( atom->getRadius() + (*it)->getRadius() ); //TODO: Speed up
		if ( dist < threshold ) {
			// If this collision is in the initial collision list, ignore it
			pair<Atom*,Atom*> collision_pair = make_pair(atom,(*it));
			if ( atom->getId() > (*it)->getId() )
				collision_pair = make_pair((*it),atom);
      auto mit=initial_collision_list.find(collision_pair);
			if ( mit!=initial_collision_list.end() ) continue;

			return true;
			//cout << "\t\t Collision: atom " << atom->getId() << " and atom " << (*it)->getId() << ", dist=" << dist << ", threshold=" << threshold << endl;
		}
	}
	return false;
}
//---------------------------------------------------------
/*
 * Given an atom and an initial list of collisions, determine if the atom is colliding with another atom
 * by ignoring atom pairs in the initial list.
 */
double Grid::minFactorWithoutCollision (Atom* atom, set< pair<Atom*,Atom*> > const &initial_collision_list, std::string collisionCheckAtoms,bool onlyCheckLargerIds) const {
	vector<Atom*> neighbors = getNeighboringAtoms(atom,onlyCheckLargerIds);
	//cout << "\t ------------------- Checking collisions for atom ... " << atom->getResidue()->getId() << " " <<  atom->getName() << " " <<  atom->Id;
	//cout << "   with " << neighbors.size() << " neighbors." << endl;
	double minFactorWithoutCollision = 999;
	for (vector<Atom*>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it) {
		if ( atom->isCovNeighbor(*it) || atom->isSecondCovNeighbor(*it) || atom->isHbondNeighbor(*it) || !(atom->isCollisionCheckAtom(collisionCheckAtoms))) {
			continue;
		}
		double dist = atom->m_Position.distanceTo((*it)->m_Position);
		double stdThreshold =  atom->getRadius() + (*it)->getRadius();
		double currentFactor = dist / stdThreshold;
		if(currentFactor < minFactorWithoutCollision ){
			// If this collision is in the initial collision list, ignore it
			pair<Atom*,Atom*> collision_pair = make_pair(atom,(*it));
			if ( atom->getId() > (*it)->getId() )
				collision_pair = make_pair((*it),atom);
//			map< pair<Atom*,Atom*>,int >::const_iterator mit=initial_collision_list.find(collision_pair);
      auto mit = initial_collision_list.find(collision_pair);
			if ( mit!=initial_collision_list.end() ) continue;

			minFactorWithoutCollision = currentFactor;
		}
	}
	return minFactorWithoutCollision;
}

//---------------------------------------------------------
/*
 * Given an atom and an initial list of collisions, get a list of colliding atoms with that atom
 * by ignoring atom pairs in the initial list.
 */
vector<Atom*> Grid::getAllCollisions (
    Atom* atom,
    set< pair<Atom*,Atom*> > const &initial_collision_list,
    std::string collisionCheckAtoms,
    bool onlyCheckLargerIds
) const
{
	vector<Atom*> collisions;
	vector<Atom*> neighbors = getNeighboringAtoms(atom,onlyCheckLargerIds);
	for (vector<Atom*>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it) {
		if ( atom->isCovNeighbor(*it) || atom->isSecondCovNeighbor(*it) || atom->isHbondNeighbor(*it) || !(atom->isCollisionCheckAtom(collisionCheckAtoms))) {
			continue;
		}
		double dist = atom->m_Position.distanceTo((*it)->m_Position);
		double threshold = m_collisionFactor * ( atom->getRadius() + (*it)->getRadius() );
		if ( dist < threshold ) {
			// If this collision is in the initial collision list, ignore it
			pair<Atom*,Atom*> collision_pair = make_pair(atom,(*it));
//			map< pair<Atom*,Atom*>,int >::const_iterator mit=initial_collision_list.find(collision_pair);
      auto mit = initial_collision_list.find(collision_pair);
			if ( mit!=initial_collision_list.end() ) continue;

			//cout << "Collision: atom1: " << atom->getResidue()->getId() << " " << atom->Name << " " << atom->Id << 
			//	      " \t atom2: " << (*it)->getResidue()->getId() << " " << (*it)->Name << " " << (*it)->Id << " , dist=" << dist << " , threshold=" << threshold << endl;
			collisions.push_back(*it);
		}
	}
	return collisions;
}
//---------------------------------------------------------

// return TRUE if the atom is removed, otherwise FALSE
bool Grid::removeAtom (Atom* atom) {
	// only delete the atom if it is indexed already
	Coordinate key = makeKey(atom->m_Position);
	map<Coordinate,vector<Atom*>*,Coordinate_cmp>::iterator map_itr=Atom_map.find(key);
	if ( map_itr!=Atom_map.end() ) {
		// search for the atom in the vector<Atom*>*
		for ( vector<Atom*>::iterator vector_itr=map_itr->second->begin(); vector_itr!=map_itr->second->end(); ++vector_itr) {
			if ( (*vector_itr)==atom ) {
				map_itr->second->erase(vector_itr);
				return true;
			}
		}
	}
	return false;
}

void Grid::addAtom (Atom* atom) {
        Coordinate key = makeKey(atom->m_Position);
        if (Atom_map.find(key)==Atom_map.end()) {
                vector<Atom*> *list = new vector<Atom*>; // this delocated? @D
                Atom_map.insert( make_pair(key,list) );
        }
        Atom_map.find(key)->second->push_back(atom);
}
