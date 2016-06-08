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
#include <stdlib.h>

#include "HbondIdentifier.h"
#include "Util.h"
#include "core/Residue.h"
#include "core/Grid.h"
#include "core/Bond.h"
#include "core/ProteinHBond.h"
#include "IO.h"
#include "math/MathUtility.h"
//#include "RunFirst.h"

using namespace std;

//RFonseca: For testing rigidity
void manualHbond(Molecule * protein, int res1, string atm1, int res2, string atm2){
	Atom *a1, *a2;
	for (vector<Atom*>::iterator atom_itr=protein->atoms.begin(); atom_itr != protein->atoms.end(); ++atom_itr) {
		Atom* a = *atom_itr;
		if(a->getResidue()->getId()==res1 && a->getName()==atm1) a1=a;
		if(a->getResidue()->getId()==res2 && a->getName()==atm2) a2=a;
	}
	Hbond * new_hb = new Hbond(a1, a2, a1->getFirstCovNeighbor(), a2->getFirstCovNeighbor());
	protein->addHbond(new_hb);
	cout<<"Manually bonded: "<<a1<<" , "<<a2<<" , "<<a1->getFirstCovNeighbor()<<" , "<<a2->getFirstCovNeighbor()<<endl;
}

void HbondIdentifier::identifyHbonds(Molecule *protein) {
	Grid grid(protein);
	for (vector<Atom*>::iterator atom_itr=protein->atoms.begin(); atom_itr != protein->atoms.end(); ++atom_itr) {
		if ( (*atom_itr)->getType() != "H" ) continue;
		// Only N (not in Pro), O, or S are allowed to be donor
		Atom* donor = (*atom_itr)->getFirstCovNeighbor();
		string donor_type = donor->getType();
		if (donor_type!="N" && donor_type!="O" && donor_type!="S")
			continue;
		if (donor_type=="N" && donor->getResidue()->getProperName()=="PRO")
			continue;
		// Check acceptor
		vector<Atom*> neighbors = grid.getNeighboringAtoms(*atom_itr,false,false,false);
		for (vector<Atom*>::iterator neighbor_itr=neighbors.begin(); neighbor_itr!=neighbors.end(); ++neighbor_itr) {
			string neighbor_type = (*neighbor_itr)->getType();
			// skip if the atom is not a possible acceptor
			if ( neighbor_type!="N" && neighbor_type!="O" && neighbor_type!="S" ) continue;
			// skip if the atom is N on main-chain
			if ( neighbor_type=="N" && !(*neighbor_itr)->isSidechainAtom()) continue;
			// skip if the atom is N TRP's side-chain
			if ( neighbor_type=="N" && (*neighbor_itr)->getResidue()->getProperName()=="TRP" && (*neighbor_itr)->isSidechainAtom() ) continue;

			// check distance H and A<2.5
			double dist_H_A = (*atom_itr)->m_Position.distanceTo( (*neighbor_itr)->m_Position);
			if ( dist_H_A>=2.5 ) continue;
			// check distance D and A<3.9
			double dist_D_A = donor->m_Position.distanceTo( (*neighbor_itr)->m_Position);
			if ( dist_D_A>=3.9 ) continue;
			// check angle D H A > 90
			double ang_D_H_A = (*atom_itr)->m_Position.getAngle( donor->m_Position, (*neighbor_itr)->m_Position);
			if ( ang_D_H_A<=90 ) continue;
			// check angle H A AA > 90
			Atom* AA = (*neighbor_itr)->getFirstCovNeighbor();
			double ang_H_A_AA = (*neighbor_itr)->m_Position.getAngle( (*atom_itr)->m_Position, AA->m_Position);
			if ( ang_H_A_AA<=90 ) continue;
			// check angle D A AA > 90
			double ang_D_A_AA = (*neighbor_itr)->m_Position.getAngle( donor->m_Position, AA->m_Position);
			if ( ang_D_A_AA<=90 ) continue;
			
			// create the Hbond object and insert into m_molecule
			//int h_id = (*atom_itr)->getId();
			//int a_id = (*neighbor_itr)->getId();
//			output << Util::i2s(h_id) << "\t" << Util::i2s(a_id) << "\t" << (*atom_itr)->getName() << " " << (*neighbor_itr)->getName() << " " << dist_H_A << " " << dist_D_A << " " << ang_D_H_A << " " << ang_H_A_AA << " " << ang_D_A_AA << endl;
//			output << Util::i2s(h_id) << "\t" << Util::i2s(a_id) << endl;
			Hbond * new_hb = new Hbond(*atom_itr, *neighbor_itr, donor, AA);
			new_hb->setIniLength(dist_H_A);
			new_hb->setIniAngle_H_A_AA(ang_H_A_AA);
			new_hb->setIniAngle_D_H_A(ang_D_H_A);
			protein->addHbond(new_hb);
			//cout<<"identifyHbonds(..) - Bonded: "<<(*atom_itr)<<" , "<<(*neighbor_itr)<<" , "<<donor<<" , "<<AA<<endl;
		}
	}
}

void HbondIdentifier::computeHbondEnergy(Molecule *protein, string path, string protein_name) {
//  //Todo: Implemeent Mayo energy function
//	IO::writeHbonds(protein,path+"hbonds.in.no_energy");
//	RunFirst::ComputeHbondEnergy(path,protein_name,"hbonds.in.no_energy","hbonds.in");
//	// get the energy values from hbonds.in
//	char temp[50];
//	int hatom_id, oatom_id;
//	double energy;
//	string atom1_sid, atom2_sid, energy_s;
//	string input_file = path + "hbonds.in";
//	ifstream input(input_file.c_str());
//	while (input.good()) {
//		input >> atom1_sid;
//		if (input.eof())
//			break;
//		input >> atom2_sid;
//		input >> energy_s;
//		input.getline(temp,50);
//		hatom_id = atoi(atom1_sid.c_str());
//		oatom_id = atoi(atom2_sid.c_str());
//		energy = atof(energy_s.c_str());
//		for (list<Hbond *>::iterator hit=protein->H_bonds.begin(); hit != protein->H_bonds.end(); ++hit) {
//			if ( (*hit)->Hatom->getId()==hatom_id && (*hit)->Acceptor->getId()==oatom_id ) {
//				(*hit)->setEnergy(energy);
//				break;
//			}
//		}
//	}
//	input.close();
}
				
void HbondIdentifier::selectHbonds(Molecule *protein, string path, string protein_name) {
	// identify hbonds if not yet
	if ( protein->H_bonds.empty() ) { identifyHbonds(protein);	}
	if ( protein->H_bonds.empty() ) { cerr<<"HbondIdentifier::selectHbonds: Was unable to identify any hydrogen bonds. Exiting."<<endl;exit(-1); }

	// compute energy if not yet
	if ( protein->H_bonds.front()->getEnergy() == DEFAULT_HBOND_ENERGY ) {
    computeHbondEnergy(protein, path, protein_name);
	}

	// select hbonds according to the decision tree
	double prob, control_prob;
	list<Hbond *>::iterator hit, hit2;
	hit = protein->H_bonds.begin();
	cout << "InfoHB)\t INITIAL H_bonds list size = " << protein->H_bonds.size() << endl;
	while ( hit!=protein->H_bonds.end() ) {
		if ( (*hit)->getLength() >= 2.4 ) {
			if ( (*hit)->getEnergy() >= 0.4 ) {
				if ( (*hit)->getLength() >= 2.67 )
					prob = 0.21;
				else
					prob = 0.37;
			}
			else {
				if ( (*hit)->getEnergy() >= -0.95 )
					prob = 0.46;
				else
					prob = 0.60;
			}
		}
		else {
			if ( (*hit)->getLength() >= 2.15 ) {
				if ( (*hit)->getLength() >= 2.27 )
					prob = 0.73;
				else
					prob = 0.84;
			}
			else {
				if ( (*hit)->getIniAngle_H_A_AA() < toRadian(104.6) )
					prob = 0.71;
				else
					prob = 0.96;
			}
		}

		prob = pow(prob,6); // 300ps

		control_prob = Random01();
		if ( prob < control_prob  ) { // need to delete this hbond
			cout << "InfoHB)\t Deleting Hbond " << (*hit)->Hatom->getId() << " " << (*hit)->Acceptor->getId() << " from hbond.in..." << endl;
			// erase from Atom's Hbond_list and Hbond_neighbor_list
			(*hit)->Hatom->removeHbond(*hit);
			(*hit)->Acceptor->removeHbond(*hit);
			// erase from Molecule's Hbnd_list
			hit2 = hit;
			++hit;
			protein->H_bonds.erase(hit2);
		}
		else {
			++hit;
		}
	}
	cout << "InfoHB)\t FINAL H_bonds list size = " << protein->H_bonds.size() << endl;
}

