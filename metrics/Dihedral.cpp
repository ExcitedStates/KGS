#include "metrics/Dihedral.h"
#include "metrics/Metric.h"
#include "core/Molecule.h"
#include "math/gsl_helpers.h"
#include <math.h>

using namespace std;

namespace metrics{

	Dihedral::Dihedral(const string& atom_selection): atom_selection(atom_selection)
	{}


	double Dihedral::distance(Configuration* c1, Configuration* c2)
	{
		int DOF = c1->getNumDOFs();

		if(DOF != c2->getNumDOFs()){
			cerr<<"Configurations to compare do not have the same number of DOFs!"<<endl;
			exit(-1);
		}
		//if(c1->getGlobalTorsions()== nullptr || c2->getGlobalTorsions()== nullptr){
		//	cerr<<"Need to set configuration and determine global torsions first!"<<endl;
		//	exit(-1);
		//}


		Molecule * m_protein = c1->getMolecule();
    bool useGlobals = m_protein!=c2->getMolecule();
    if(useGlobals){
      c1->updateGlobalTorsions();
      c2->updateGlobalTorsions();
    }


		int count = 0;
		double distance=0.0;
    for(KinEdge*& edge: m_protein->m_spanning_tree->Edges){
			int dofId = edge->getDOF()->getIndex();
			double angle_diff;
			if(useGlobals) {
        angle_diff = fabs(c2->getGlobalTorsion(dofId) - c1->getGlobalTorsion(dofId));
        if(angle_diff>M_PI)
          angle_diff = 2*M_PI - angle_diff;
      }else{
        angle_diff = fabs(c2->m_dofs[dofId] - c1->m_dofs[dofId]);
        if(angle_diff>M_PI)
          angle_diff = 2*M_PI - angle_diff;
      }
      //cout<<"DOFid: "<<dofId<<" , diff: "<<angle_diff<<endl;

      //if(edge->getBond()) {
      //  Atom *a=edge->getBond()->Atom1;
      //  Atom *b=edge->getBond()->Atom2;
      //  if( a->getResidue()->getId() == 16 && a->getName() == "CA" ) {
      //    cout << "Dihedral::Distance - " << *a <<" - "<<*b<< " : " << angle_diff << endl;
      //  }
      //}

//      cout<<"Dihedral::Distance - "<< c2->getGlobalTorsion(dofId)<<"-"<< c1->getGlobalTorsion(dofId)<<endl;
			if(atom_selection=="MOV"){
				int cycle_dof_id = edge->getDOF()->getCycleIndex();
				bool locked = edge->getBond()->constrained;
				if( cycle_dof_id == -1 ){//free dihedral, always moveable
					distance += angle_diff*angle_diff;
//          if(angle_diff>0.0001)
            count++;
				}
				else if (!locked ) {
					distance += angle_diff*angle_diff;
//          if(angle_diff>0.0001)
					count++;
				}
			}
			else{
				distance += angle_diff*angle_diff;
//        if(angle_diff>0.0001)
				count++;
			}
//			if( Abs(angle_diff - diff_rel ) > 0.0001)
//				cout<<"Difference at "<<dofId<<": "<<angle_diff - diff_rel<<" between atoms "<<(*eit)->getBond()->Atom1<<" and "<<(*eit)->getBond()->Atom2<<endl;;
		}

//		cout<<"relative distance: "<<sqrt(distanceRel)<<" absolute distance: "<<sqrt(distance)<<endl;
//		return sqrt(distance/count);
		return sqrt(distance/count);
	}
}
