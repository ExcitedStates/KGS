#include "metrics/Dihedral.h"
#include "metrics/Metric.h"
#include "core/Molecule.h"

using namespace std;

namespace metrics{

	Dihedral::Dihedral(const string& atom_selection): atom_selection(atom_selection)
	{}


	double Dihedral::distance(Configuration* c1, Configuration* c2)
	{
		int DOF = c1->getNumDOFs();

		if(DOF != c2->getNumDOFs()){
			cerr<<"Configurations to compare do not have the same number of m_dofs!"<<endl;
			exit(-1);
		}
		//if(c1->getGlobalTorsions()== nullptr || c2->getGlobalTorsions()== nullptr){
		//	cerr<<"Need to set configuration and determine global torsions first!"<<endl;
		//	exit(-1);
		//}

		Molecule * m_protein = c1->getMolecule();

		int count = 0;
		double distance=0.0, distanceRel=0.0;
    for(KinEdge*& edge: m_protein->m_spanning_tree->Edges){
			int dofId = edge->getDOF()->getIndex();
			//double angle_diff = c2->getGlobalTorsions(dofId)-c1->getGlobalTorsions(dofId);
			//angle_diff = formatRangeRadian(angle_diff);
			double diff_rel = formatRangeRadian(c2->m_dofs[dofId] - c1->m_dofs[dofId]);
			if(atom_selection=="MOV"){
				int cycle_dof_id = edge->getDOF()->getCycleIndex();
				bool locked = edge->getBond()->constrained;
				if( cycle_dof_id == -1 ){//free dihedral, always moveable
					//distance += angle_diff*angle_diff;
					distanceRel += diff_rel*diff_rel;
					count++;
				}
				else if (!locked ) {
					//distance += angle_diff*angle_diff;
					distanceRel += diff_rel*diff_rel;
					count++;
				}
			}
			else{
				//distance += angle_diff*angle_diff;
				distanceRel += diff_rel*diff_rel;
				count++;
			}
//			if( Abs(angle_diff - diff_rel ) > 0.0001)
//				cout<<"Difference at "<<dofId<<": "<<angle_diff - diff_rel<<" between atoms "<<(*eit)->getBond()->Atom1<<" and "<<(*eit)->getBond()->Atom2<<endl;;
		}

//		cout<<"relative distance: "<<sqrt(distanceRel)<<" absolute distance: "<<sqrt(distance)<<endl;
//		return sqrt(distance/count);
		return sqrt(distanceRel/count);
	}
}
