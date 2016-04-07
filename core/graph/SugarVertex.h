

#ifndef SUGARVERTEX_H
#define SUGARVERTEX_H

#include <iostream>
#include <map>

#include "KinGraph.h"
#include "math/MathUtility.h"
#include "core/Bond.h"
#include "core/Rigidbody.h"
#include "math3d/primitives.h"

using namespace Math3D;


class SugarVertex: public KinVertex{
	public:
		int DOF_id;//Can be different from -1 if its a sugar
		int Cycle_DOF_id;

		/* Array of all atoms in this graph vertex. Each entry holds one sugar-ring atom member (first list-element) and atoms not attached
		 * to any other member. The order is such that adjacent sugar-ring atom members are adjacent in the array and the vector from member
		 * 0 to member 1 is the main torsion angle. */
		std::list<Atom*> ringAttached[5];
		double angles[5];//Angles between ring-members 4-0-1, 0-1-2, 1-2-3, 2-3-4 and 3-4-0
		double dists[5];//Distances between ring-members 0-1, 1-2, 2-3, 3-4, 4-0
		double initTorsion;//Initial torsion angle of bond 0-1
		double initUp;
		double Amplitude;//Largest torsion of bond 0-1 such that the sugar can still close.

		SugarVertex ();
		SugarVertex (int id, Rigidbody* rb);
		~SugarVertex ();

		void setParent(KinVertex* v);

		double getCurrentUP();
		void configurationToGlobalMatrix(RigidTransform* ms, double* m_f);//, bool usePosition2);
		void updateDelayed();//bool usePosition2);
		Math3D::Vector3 computeJacobianEntry(Edge* outEdge, double* m_f, Vector3& p);
	private:
		bool getSugarConformation(double tau, Vector3* conf, double Am);//, bool usePosition2 = false);
		Edge* getEdgeStartingAt(std::string atomName);
		void collectAttached(Atom* a, Rigidbody* rb, std::list<Atom*>& visited);
		void collectRest(Rigidbody* rb);
		bool sugarMember(Atom* a);
		std::map<int,Math3D::Vector3> delayedUpdates;
};

#endif
