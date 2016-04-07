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
#include "string.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <new>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <metrics/RMSD.h>


#include "CTKTimer.h"
#include "IO.h"
#include "core/Atom.h"
#include "core/Chain.h"
#include "core/Bond.h"
#include "core/ProteinHBond.h"
#include "DisjointSets.h"
#include "Logger.h"
#include "Selection.h"
#include "Color.h"


using namespace std;

ResidueProfile IO::readResidueProfile () {

	//Read residue profile from the array in ResidueProfiles.h
	ResidueProfile residue_profile;
	map< string,vector<CovBond> >::iterator it;

	for(int rp=0;;rp++){
		if(strcmp(COV_BOND_PROFILES[rp][0],"END")==0) break;

		it = residue_profile.find(COV_BOND_PROFILES[rp][0]);
		if(it==residue_profile.end()) {
			vector<CovBond> cov_bond_v;
			it = residue_profile.insert(residue_profile.begin(), make_pair(COV_BOND_PROFILES[rp][0], cov_bond_v));
		}
		
		CovBond cov_bond(COV_BOND_PROFILES[rp][1], COV_BOND_PROFILES[rp][2]);
		it->second.push_back(cov_bond);
	}

	return residue_profile;
}

void IO::readPdb (Molecule * protein, string pdb_file, vector<string> &extraCovBonds, Molecule * reference) {
	// Read all atoms. Create one atom for each ATOM record in file
	if ( protein->getName().compare("UNKNOWN")==0 ) {
		protein->setName(pdb_file);
	}
	ifstream pdb(pdb_file.c_str());
	if (!pdb.good()) {
		cerr << "Error: cannot open file " << pdb_file << endl;
		exit(1);
	}
	string line, temp;
	while (!pdb.eof()) {
		getline(pdb,line);
		if (pdb.eof()) break;
		if (line.substr(0,4)!="ATOM") // skip if it is not an ATOM line
			continue;
		// chain info
		string chain_name = line.substr(21,1); // line[22]
		// residue info
		int res_id = atoi(line.substr(22,4).c_str()); // line[23:26]
		string res_name = Util::trim(line.substr(17,3)); // line[18:20]
		// atom info
		int atom_id = atoi(line.substr(6,5).c_str()); // line[7:11]s
		string atom_name = Util::trim(line.substr(12,5)); // line[13:17]
		if(atom_name=="OP3") continue;
		if (atom_name.at(0)>=49 && atom_name.at(0)<=57) { // if the first char is 1-9
			string temp_name(atom_name.substr(1,3));
			temp_name += atom_name.substr(0,1);
			atom_name = temp_name;
		}

		double x = atof(line.substr(30,8).c_str()); // line[31:38]
		double y = atof(line.substr(38,8).c_str()); // line[39:46]
		double z = atof(line.substr(46,8).c_str()); // line[47:54]

		Coordinate pos(x,y,z);
		protein->addAtom(chain_name,res_name,res_id,atom_name,atom_id,pos);
	}
	pdb.close();

	ResidueProfile residue_profile = readResidueProfile();

	for (auto const& cur_chain: protein->chains) {
		for (auto const& cur_res: cur_chain->getResidues()){
			string res_name = cur_res->getName();

			map< string, vector<CovBond> >::iterator profile_it=residue_profile.find(res_name);
			if (profile_it==residue_profile.end()) {
				cerr << "No such residue (" << res_name << ") in profile." << endl;
				exit(1);
			}
			vector<CovBond> bonds = profile_it->second;
			for (vector<CovBond>::iterator bond_it=bonds.begin(); bond_it!=bonds.end(); ++bond_it) {
				string atom_name1 = bond_it->first;
				string atom_name2 = bond_it->second;
				Residue* res1 = cur_res;
				Residue* res2 = cur_res;
				if (atom_name1[0]=='-') {
					atom_name1 = Util::trim(atom_name1,'-');
					res1 = cur_res->getLastResidue();
				}
				else if (atom_name1.at(0)=='+') {
					atom_name1 = Util::trim(atom_name1,'+');
					res1 = cur_res->getNextResidue();
				}
				if (atom_name2.at(0)=='-') {
					atom_name2 = Util::trim(atom_name2,'-');
					res2 = cur_res->getLastResidue();
				}
				else if (atom_name2.at(0)=='+') {
					atom_name2 = Util::trim(atom_name2,'+');
					res2 = cur_res->getNextResidue();
				}
				if (res1==NULL || res2==NULL) {
					continue;
				}

				makeCovBond(res1,res2,atom_name1,atom_name2);
			} // finish looping over bonds
		} // finish looping over residues
	} // finish looping over chains

	if(reference == NULL) {
		for (unsigned int i = 0; i < extraCovBonds.size(); i++) {
			vector<string> tokens = Selection::split(extraCovBonds[i], "-");
			int atomId1 = atoi(tokens[0].c_str());
			int atomId2 = atoi(tokens[1].c_str());
			Atom *a1 = protein->getAtom(atomId1);
			Atom *a2 = protein->getAtom(atomId2);
			if (a1 == NULL) {
				cerr << "Cannot find atom with id " << atomId1 << endl;
				exit(-1);
			}
			if (a2 == NULL) {
				cerr << "Cannot find atom with id " << atomId2 << endl;
				exit(-1);
			}
			makeCovBond(a1->getResidue(), a2->getResidue(), a1->getName(), a2->getName());
			cout << "Creating bond between " << a1 << " and " << a2 << " in protein " << protein->getName() << endl;
		}
	}
	else {
		for (unsigned int i = 0; i < extraCovBonds.size(); i++) {
			vector<string> tokens = Selection::split(extraCovBonds[i], "-");
			int atomId1 = atoi(tokens[0].c_str());
			int atomId2 = atoi(tokens[1].c_str());
			Atom* a1 = reference->getAtom(atomId1);//id's from reference
			Atom* a2 = reference->getAtom(atomId2);
			//residue Ids and names
			int resId1 = a1->getResidue()->getId();
			string name1 = a1->getName();
			string chainName1 = a1->getResidue()->getChain()->getName();
			int resId2 = a2->getResidue()->getId();
			string name2 = a2->getName();
			string chainName2 = a2->getResidue()->getChain()->getName();
			//use names to identify atoms in protein
			Atom* a3 = protein->getAtom(chainName1,resId1, name1);
			Atom* a4 = protein->getAtom(chainName2,resId2, name2);
			if(a3==NULL) { cerr<<"Cannot find atom with residue id "<<resId1<<" and name "<<name1<<endl; exit(-1); }
			if(a4==NULL) { cerr<<"Cannot find atom with residue id "<<resId2<<" and name "<<name2<<endl; exit(-1); }
			makeCovBond(a3->getResidue(), a4->getResidue(), a3->getName(), a4->getName());		}
	}

	// Fill in the second_cov_neighbor_list
	// Cannot do this step when the bond is still in creation because it won't know the neighbors of its neighbors yet
	for (vector<Atom*>::iterator ait=protein->atoms.begin(); ait != protein->atoms.end(); ++ait) {
		for (vector<Atom*>::iterator n1=(*ait)->Cov_neighbor_list.begin(); n1!=(*ait)->Cov_neighbor_list.end(); ++n1) {
			for (vector<Atom*>::iterator n2=(*n1)->Cov_neighbor_list.begin(); n2!=(*n1)->Cov_neighbor_list.end(); ++n2) {
				// check if n2 is ait itself. if yes, ignore it.
				if ( (*n2) == (*ait) )
					continue;
				// check whether n2 is already in the second_cov_neighbor_list of ait
				bool got_already = false;
				for (vector<Atom*>::iterator n3=(*ait)->Second_cov_neighbor_list.begin(); n3!=(*ait)->Second_cov_neighbor_list.end(); ++n3) {
					if ( (*n3) == (*n2) ) {
						got_already = true;
						break;
					}
				}
				if (!got_already)
					(*ait)->Second_cov_neighbor_list.push_back(*n2);
			}
		}

		//Print warning if atom has no covalent neighbors
		if( (*ait)->Cov_neighbor_list.size()==0 ){
			cerr<<"IO::readPdb - Error: Atom "<<(*ait)<<" has no covalent neighbors. Probably the atom-name isn't in the residue profiles"<<endl;
			exit(-1);
		}
	}

	// Change the number of bars of locked bonds and peptide bonds to 6.
	// Cannot do this step when the bond is just created because it won't know the number of neighbors yet
	// No need to work on H-bonds because we assume all of them have 5 bars.
	for (list<Bond *>::iterator bit=protein->Cov_bonds.begin(); bit != protein->Cov_bonds.end(); ++bit) {
		if ( (*bit)->isLocked() || (*bit)->Bond::isPeptideBond() )
			(*bit)->Bars = 6;
	}

	// Assign main-chain/side-chain for all atoms
	for (vector<Atom*>::iterator it=protein->atoms.begin(); it != protein->atoms.end(); ++it) {
		string atom_name = (*it)->getName();
		if (atom_name.compare("N")==0 || atom_name.compare("CA")==0 || atom_name.compare("C")==0 || atom_name.compare("O")==0) // this is a main-chain atom
			(*it)->setAsMainchainAtom();
		else if (
				atom_name.compare("P")==0 || atom_name.compare("OP1")==0 || atom_name.compare("OP2")==0 ||
				atom_name.compare("O5'")==0 ||
				atom_name.compare("C5'")==0 ||
				atom_name.compare("C4'")==0 || atom_name.compare("O4'")==0 ||
				atom_name.compare("C3'")==0 ||
				atom_name.compare("O3'")==0 ) // this is a main-chain atom in RNA
			(*it)->setAsMainchainAtom();
		else if (atom_name.find("H")!=string::npos) { // this is a H-atom
			if((*it)->Cov_neighbor_list.size() > 0) {
				Atom* cov = (*it)->getIthCovNeighbor(0);
				string cov_name = cov->getName();
				if (cov_name.compare("N")==0 || cov_name.compare("CA")==0 || cov_name.compare("C")==0 || cov_name.compare("O")==0 ||
						cov_name.compare("P")==0 || cov_name.compare("OP1")==0 || cov_name.compare("OP2")==0 ||
						cov_name.compare("O5'")==0 ||
						cov_name.compare("C5'")==0 ||
						cov_name.compare("C4'")==0 ||	cov_name.compare("O4'")==0 ||
						cov_name.compare("C3'")==0 ||
						cov_name.compare("O3'")==0 )
					(*it)->setAsMainchainAtom();
			}
			else
				(*it)->setAsSidechainAtom();
		}
		else
			(*it)->setAsSidechainAtom();
	}

	// Index m_protein atoms
	protein->indexAtoms();
	protein->backupAtomIndex();

//		m_protein->printSummaryInfo();
}

void IO::makeCovBond (Residue* res1, Residue* res2, string atom_name1, string atom_name2) {
  //TODO: None of this needs to be in IO. Create a function in Molecule that does this and call it directly
	Atom* atom1 = res1->getAtom(atom_name1);
	Atom* atom2 = res2->getAtom(atom_name2);
	bool valid = true;
	if (atom1==NULL) {
		//cerr << "Warning: No such atom (" << atom_name1 << ") in Res" << res1->getId() << endl;
		valid = false;
	}
	if (atom2==NULL) {
		//cerr << "Warning: No such atom (" << atom_name2 << ") in Res" << res2->getId() << endl;
		valid = false;
	}
	if (valid) {
		Bond * new_cb = new Bond(atom1, atom2, "COV");
		res1->getChain()->getProtein()->addCovBond(new_cb);
	}
}

void IO::readDssp (Molecule * protein, string dssp_file) {
	ifstream input(dssp_file.c_str());
	if (!input.good()) {
		cerr << "Error: cannot open file " << dssp_file << endl;
		exit(1);
	}
	string sse_symbol="", last_symbol="x", chain_name = "", last_chain_name="";
	int res_id, last_res_id=-100;
	Chain* chain;
	Residue* res;
	string line;
	bool start_read = false;
	int helix_index=0, loop_index=0;
	while ( getline(input,line) ) {
		if (input.eof()) break;
		if (line.at(2)=='#') {start_read = true; continue;}
		if (!start_read) continue;
		if (line.at(13)=='!') continue; // there is a discontinuity in either the backbone or the chain. skip this line

		// Read in the symbol in DSSP file
		chain_name = line.substr(11,1); // line[12]
		if (chain_name.compare(last_chain_name)!=0) {
			chain = protein->getChain(chain_name);
			if (chain==NULL) {
				cerr << "Error: No such chain (" << chain_name << ") in m_protein." << endl;
				exit(1);
			}
			last_chain_name = chain_name;
		}
		res_id = atoi(line.substr(6,5).c_str()); // line[7:11]
		res = chain->getResidue(res_id);
		if (res==NULL) {
			cerr << "Error: No such residue " << res_id << " in chain (" << chain_name << ")." << endl;
			exit(1);
		}
		sse_symbol = line.substr(16,1); // line[17]
		if (sse_symbol.compare(" ")==0) // represents for loop, change the symbol to L
			sse_symbol = "L";

		// Convert the symbols and number the helices and loops
		// Currently, no numbering on strands because no information on how to pair strands into sheets.
		// G, H, I are all helices -> A. 
		// E and B are strands -> B.
		// The rest are loops -> L.
		if (sse_symbol.compare("G")==0 || sse_symbol.compare("H")==0 || sse_symbol.compare("I")==0) { // this is a helix
			if (last_symbol.at(0)!='A' || last_res_id!=res_id-1)
				++helix_index;
			sse_symbol = "A"+Util::i2s(helix_index);
		}
		else if (sse_symbol.compare("E")==0 || sse_symbol.compare("B")==0) { // this is a strand
			sse_symbol = "B";
		}
		else {
			if (last_symbol.at(0)!='L' || last_res_id!=res_id-1)
				++loop_index;
			sse_symbol = "L"+Util::i2s(loop_index);
		}

		// Assign to the residue
		res->setSSE(sse_symbol);
		last_res_id = res_id;
		last_symbol = sse_symbol;
	}
	input.close();
}


void IO::readRigidbody (Molecule * protein) {
	//Create disjoint set
	DisjointSets ds(protein->atoms[protein->size() - 1]->getId() + 1); //Assumes the last atom has the highest id.

	//For each atom, a1, with exactly one cov neighbor and not participating in an hbond, a2, call Union(a1,a2)
	for (int i=0;i<protein->size();i++){
		Atom* atom = protein->atoms[i];
		if(atom->Cov_neighbor_list.size()==1 && atom->Hbond_neighbor_list.size()==0){
			ds.Union(atom->getId(), atom->Cov_neighbor_list[0]->getId());
			//cout<<"Only one neighbor: "<<atom->getName()<<" "<<atom->getId()<<" - "<<atom->Cov_neighbor_list[0]->getName()<<" "<<atom->Cov_neighbor_list[0]->getId()<<endl;
		}
	}

	//Create a map to determine if a covalent bond is fixed
	set<string> fixedBonds;
	for(int rp=0;;rp++){
		if(strcmp(FIXED_BOND_PROFILES[rp][0],"END")==0) break;
		string tmp(FIXED_BOND_PROFILES[rp][0]);
		tmp+="_";
		string atomName1(FIXED_BOND_PROFILES[rp][1]);
		string atomName2(FIXED_BOND_PROFILES[rp][2]);
		if(atomName1<atomName2)	tmp+=atomName1+"_"+atomName2;
		else					tmp+=atomName2+"_"+atomName1;
        //log("debugRas")<<"Inserting "<<tmp<<endl;;
		fixedBonds.insert(tmp);
	}

	//For each fixed bond (a1,a2) call Union(a1,a2)
	for (list<Bond *>::iterator it=protein->Cov_bonds.begin(); it != protein->Cov_bonds.end(); ++it){
		Bond * bond = *it;

		///****************************************
		///UPDATED: we use two ways for now to identify locked bonds!
		///The isLocked/isPeptide bond check in Bond and the profiles
		///As for some profiles, it can be either locked or rotatable!

		if( bond->Bars == 6){//This is fixed in the Bond -> isLocked
      //cout<<"IO::readRigidBody - Bars=6"<<endl;
			ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
			continue;
		}
		
		Atom* a1 = bond->Atom1;
		Atom* a2 = bond->Atom2;
		//If the residue containing a1 (R1) precedes the one containing a2 (R2) there are four patterns in 
		//FIXED_BOND_PROFILES that can match: R1_a1_+a2, R1_+a2_a1, R2_-a1_a2 and R2_a2_-a1. Since atom names
		//were ordered in the previous loop, however, we only have to check R1_+a2_a1 and R2_-a1_a2. 
		string query1, query2;
		if(  a1->getResidue()->getNextResidue() == a2->getResidue()  ){
			query1 = a1->getResidue()->getName()+"_+"+a2->getName()+"_"+a1->getName();
			query2 = a2->getResidue()->getName()+"_-"+a1->getName()+"_"+a2->getName();
		}else if(  a2->getResidue()->getNextResidue() == a1->getResidue()  ){
			query1 = a1->getResidue()->getName()+"_+"+a1->getName()+"_"+a2->getName();
			query2 = a2->getResidue()->getName()+"_-"+a2->getName()+"_"+a1->getName();
		}else{
			//If a1 and a2 are in the same residue simply put the lexicographically smallest name first. 
			if(a1->getName()<a2->getName())
				query1 = a1->getResidue()->getName()+"_"+a1->getName()+"_"+a2->getName();
			else
				query1 = a1->getResidue()->getName()+"_"+a2->getName()+"_"+a1->getName();
		}
        //log("debugRas")<<bond<<" .. query: "<<query1<<" ("<<query2<<")"<<endl;

		if(fixedBonds.find(query1)!=fixedBonds.end() || fixedBonds.find(query2)!=fixedBonds.end()){
			ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
			///Set the number of bars to six for all these bonds!
			(*it)->Bars = 6;
		}
	}


	int c=0;
	map<int,int> idMap;//Maps atom id's to rigid body id's for use in the DS structure.

	//Map the set-ID's to RB-ID's and add bonded atoms to RBs.
	for (int i=0;i<protein->size();i++){
		Atom* atom = protein->atoms[i];

		//Map the set-id to the RB-id
		int set_id = ds.FindSet(atom->getId());
		int body_id;
		if(idMap.find(set_id)!=idMap.end()) body_id = idMap.find(set_id)->second;
		else {
			body_id = c++;
			idMap.insert( make_pair(set_id, body_id) );
		}
		//If the set containing a1 is not a rigid body: create one
		if ( protein->Rigidbody_map_by_id.find(body_id)==protein->Rigidbody_map_by_id.end() ) {
			Rigidbody* new_rb = new Rigidbody(body_id);
			protein->Rigidbody_map_by_id.insert( make_pair(body_id,new_rb) );
		}
		Rigidbody* rb = protein->Rigidbody_map_by_id.find(body_id)->second;
		if (!rb->containsAtom(atom)) rb->addAtom(atom);

		// Add the covalent neighbors of this atom into its rigid body, if it's not there already
//		for ( vector<Atom*>::iterator it=atom->Cov_neighbor_list.begin(); it!=atom->Cov_neighbor_list.end(); ++it) {
//			if ( !rb->containsAtom(*it) ) {
//				rb->addAtom(*it);
//			}
//		}
		//Add the Hbond neighbors of this atom into this rigid body, if it's not there
//		for ( vector<Atom*>::iterator it=atom->Hbond_neighbor_list.begin(); it!=atom->Hbond_neighbor_list.end(); ++it) {
//			if ( !rb->containsAtom(*it) ) {
//				rb->addAtom(*it);
//			}
//		}
	}

	//for(map<unsigned int, Rigidbody*>::iterator it = m_protein->Rigidbody_map_by_id.begin(); it!=m_protein->Rigidbody_map_by_id.end(); it++){
	//	Rigidbody* rb = it->second;
    //	log("debugRas")<<rb<<endl;
	//}

	//Delete small RBs and sort atoms within each RB
	map<unsigned int, Rigidbody*>::iterator it = protein->Rigidbody_map_by_id.begin();
	while ( it!=protein->Rigidbody_map_by_id.end() ) {
		Rigidbody *rb = it->second;
		if ( rb->Atoms.size() <= 0 ) {
			cerr << "Error: rigid body " << rb->id() << " has less than 3 atoms." << endl;
			//exit(1);
			//for (vector<Atom*>::iterator ait=rb->Atoms.begin(); ait!=rb->Atoms.end(); ++ait) {
      //  (*ait)->removeRigidbody(rb);
			//}
			protein->Rigidbody_map_by_id.erase(it++);
			delete rb;
		}
		else { // sort atom ids
#ifndef WIN32 // Liangjun Zhang's tmp code
			vector<Atom*>::iterator sit = rb->Atoms.begin();
			vector<Atom*>::iterator eit = rb->Atoms.end();

			sort(sit,eit,Atom::compare);
#endif
			// Determine if Rigidbody is on Mainchain
			rb->setMainchainRb();
			++it;
		}
	}
	
	//Store bonds in rigid bodies
	for (list<Bond *>::iterator bit=protein->Cov_bonds.begin(); bit != protein->Cov_bonds.end(); ++bit) {
		int setId1 = ds.FindSet((*bit)->Atom1->getId());
		protein->Rigidbody_map_by_id.find( idMap.find(setId1)->second )->second->addBond(*bit);
		int setId2 = ds.FindSet((*bit)->Atom2->getId());
		if(setId1!=setId2)
			protein->Rigidbody_map_by_id.find( idMap.find(setId2)->second )->second->addBond(*bit);
	}
	for (list<Hbond *>::iterator bit=protein->H_bonds.begin(); bit != protein->H_bonds.end(); ++bit) {
		int setId1 = ds.FindSet((*bit)->Atom1->getId());
		protein->Rigidbody_map_by_id.find( idMap.find(setId1)->second )->second->addBond(*bit);
		int setId2 = ds.FindSet((*bit)->Atom2->getId());
		if(setId1!=setId2)
			protein->Rigidbody_map_by_id.find( idMap.find(setId2)->second )->second->addBond(*bit);
	}


}

void IO::writePdb (Molecule * protein, string output_file_name) {
	for(int i=0;i<output_file_name.length();i++){
		if(output_file_name[i]=='/'){
			mkdir(output_file_name.substr(0,i).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
	}
	ofstream output(output_file_name.c_str());
	if(!output.is_open()) {
		//cerr<<"Cannot write to "<<output_file_name<<". You might need to create output directory first"<<endl;
		cerr<<"Cannot write to "<<output_file_name<<endl;
		exit(-1);
	}
	if(protein->m_conf!=NULL){
		Configuration* c = protein->m_conf;
		output << "REMARK\tID = " << c->m_id << endl;
		if(c->getParent() != NULL)
			output << "REMARK\tParent ID = " << c->getParent()->m_id << endl;
		output << "REMARK\tTree depth = " << c->m_treeDepth << endl;
		output<<"REMARK\tTree-path = ";
		int count=0;
		while(c->getParent() != NULL && c->getParent() != c) { output << c->m_id << " "; c=c->getParent(); count++; }
		output << c->m_id << endl;
		output<<"REMARK\tDistance_initial = "<<setprecision(3)<<protein->m_conf->m_distanceToIni<<endl;
		output<<"REMARK\tDistance to m_parent = "<<setprecision(3)<<protein->m_conf->m_distanceToParent<<endl;
		output<<"REMARK\tMax violation = "<<setprecision(3)<<protein->m_conf->m_maxConstraintViolation<<endl;
	}
	for (vector<Atom*>::iterator atom_itr=protein->atoms.begin(); atom_itr != protein->atoms.end(); ++atom_itr) {
		Atom* atom = *atom_itr;
		Residue* res = (*atom_itr)->getResidue();
		char buffer[100];
		sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  ",
		//sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6d          %2s  ",
			atom->getId(),atom->getName().c_str(),
			res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
			atom->m_Position.x,atom->m_Position.y,atom->m_Position.z,atom->getType().c_str());
			//atom->m_Position.x,atom->m_Position.y,atom->m_Position.z,atom->getBiggerRigidbody()->id(),atom->getType().c_str() );
		string line(buffer);
		output << line << endl;
	}
	output.close();
}

void IO::writeQ (Molecule *protein, Configuration* referenceConf, string output_file_name) {
	ofstream output(output_file_name.c_str());
	if(!output.is_open()) {
		cerr<<"Cannot write to "<<output_file_name<<endl;
		exit(-1);
	}
	for (int i=0; i<protein->totalDofNum(); i++) {
		output << protein->m_conf->m_dofs[i] << "  " << protein->m_conf->getGlobalTorsions(i) << "  " << protein->m_conf->m_sumProjSteps[i] << endl;
//		output <<  m_protein->m_conf->m_f[i] << "  " << m_protein->m_conf->getGlobalTorsions(i)<< endl;
	}
	output.close();

	///For testing only
	vector<Edge*>::iterator eit;
	int i=0;
	ofstream myfile;
	string myfile_s = output_file_name.substr(0,output_file_name.length()-4) + "_absChange.txt";
	myfile.open(myfile_s.c_str());
	//myfile<<"Absolute change, dof_id, res_id, 0/1/2 for free/coordinated/rigid"<<endl;
	//Keep track of changes (in magnitude)
	for (eit = protein->m_spanning_tree->Edges.begin(); eit != protein->m_spanning_tree->Edges.end(); eit++){
		int dof_id = (*eit)->DOF_id;
		int resId = (*eit)->getBond()->Atom1->getResidue()->getId();
//		if( dof_id == -1)
//			continue;
//		else{
			int cycleDOF_id = (*eit)->Cycle_DOF_id;
			int second = 1;
			if(cycleDOF_id == -1)
				second = 0;
			//else if(gsl_vector_get(m_protein->m_conf->CycleNullSpace->m_rigidAngles,cycleDOF_id)==1)
      else if( protein->m_conf->getNullspace()->IsAngleRigid(cycleDOF_id) )
				second = 2;
			double absChangeI = formatRangeRadian(protein->m_conf->getGlobalTorsions(i) - referenceConf->getGlobalTorsions(i));
			myfile << absChangeI <<" "<< dof_id <<" "<<resId <<" "<<second<<endl;
			i++;
//		}
	}
	myfile.close();

}

void IO::writeBondLengthsAndAngles (Molecule *molecule, string output_file_name) {
	string covBondFile = output_file_name + "_allCovBonds.txt";
	ofstream output(covBondFile.c_str());
	for (list<Bond*>::iterator bond_itr=molecule->Cov_bonds.begin(); bond_itr!=molecule->Cov_bonds.end(); ++bond_itr) {
		Bond* bond = (*bond_itr);
		Vector3 bondVec = bond->Atom1->m_Position - bond->Atom2->m_Position;
		output << bondVec.norm() << endl;
	}
	output.close();

	string edgeFile = output_file_name + "_allEdges.txt";
	ofstream output1(edgeFile.c_str());
	for (vector<Edge*>::iterator edge_itr=molecule->m_spanning_tree->Edges.begin(); edge_itr!=molecule->m_spanning_tree->Edges.end(); ++edge_itr) {
		Bond* bond = (*edge_itr)->getBond();
		Vector3 bondVec = bond->Atom1->m_Position - bond->Atom2->m_Position;
		output1 << bondVec.norm() << endl;
	}
	output1.close();

	string anchorEdgeFile = output_file_name + "_allAnchors.txt";
	ofstream output2(anchorEdgeFile.c_str());
	for (vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator edge_itr=molecule->m_spanning_tree->CycleAnchorEdges.begin(); edge_itr!=molecule->m_spanning_tree->CycleAnchorEdges.end(); ++edge_itr) {
		Bond* bond = edge_itr->first->getBond();
		Vector3 bondVec = bond->Atom1->m_Position - bond->Atom2->m_Position;
		output2 << bondVec.norm() << endl;
	}
	output2.close();
}


void IO::writeCovBonds (Molecule *protein, string output_file_name) {
	ofstream output(output_file_name.c_str());
	for (list<Bond *>::iterator bond_itr=protein->Cov_bonds.begin(); bond_itr != protein->Cov_bonds.end(); ++bond_itr) {
		output << right << setw(8) << (*bond_itr)->Atom1->getId() << right << setw(8) << (*bond_itr)->Atom2->getId() << right << setw(8) << (*bond_itr)->Bars << endl;
	}
	output.close();
}

void IO::readCovBonds (Molecule *protein, string in_file_name) {

	ifstream cov(in_file_name.c_str());
	if (!cov.good()) {
		cerr << "Error: cannot open file " << in_file_name << endl;
		exit(1);
	}
	string line, temp;
	while (!cov.eof()) 
	{
		getline(cov,line);
        // log("debugRas") << line;
		if(line.size()== 0) // the last empty line
			break;

		string atom_name1 = line.substr(0,8);
		string atom_name2 = line.substr(8,8);
		Atom *atom1 = protein->getAtom(atoi(atom_name1.c_str()));
		Atom *atom2 = protein->getAtom(atoi(atom_name2.c_str()));

		Bond * new_cb = new Bond(atom1, atom2, "COV");
		protein->addCovBond(new_cb);
	}
	cov.close();
}

void IO::writeHbonds (Molecule *protein, string output_file_name) {
	ofstream output(output_file_name.c_str());
	for (list<Hbond *>::iterator hb_itr=protein->H_bonds.begin(); hb_itr != protein->H_bonds.end(); ++hb_itr) {
		output << std::right << setw(8) << (*hb_itr)->Hatom->getId();
		output << std::right << setw(8) << (*hb_itr)->Acceptor->getId();
		if ( (*hb_itr)->Energy==DEFAULT_HBOND_ENERGY ) {
			output << endl;
		}
		else {
			output << std::right << setw(16) << (*hb_itr)->Energy;
			output << std::right << setw(5) << 5;
			output << std::right << setw(16) << (*hb_itr)->Dist_H_A;
			output << std::right << setw(16) << (*hb_itr)->Ang_H_A_AA << endl;
		}
	}
	output.close();
}

void IO::readHbonds (Molecule *protein, string hbond_file_name) {
	ifstream input(hbond_file_name.c_str());
	std::string line;
	string atom1_sid, atom2_sid, energy_s;
	int hatom_id, oatom_id;
	double energy;
	Atom *hatom, *oatom, *donor, *AA;
	while(std::getline(input,line)){
		std::istringstream input(line);
		input >> atom1_sid;
		input >> atom2_sid;
		if (input >> energy_s){
			energy = atof(energy_s.c_str());
		}else{
			energy = DEFAULT_HBOND_ENERGY;
		}
		hatom_id = atoi(atom1_sid.c_str());
		oatom_id = atoi(atom2_sid.c_str());

		hatom = protein->getAtom(hatom_id);
        if(hatom==NULL){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<hatom_id<<endl; exit(-1); }
        oatom = protein->getAtom(oatom_id);
        if(oatom==NULL){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<oatom_id<<endl; exit(-1); }

		donor = hatom->getFirstCovNeighbor();
		AA = oatom->getFirstCovNeighbor();
		Hbond * new_hb = new Hbond(hatom, oatom, donor, AA, energy);
		protein->addHbond(new_hb);
	}
	input.close();
}

void IO::readAnnotations (Molecule *protein, string annotation_file_name){
    ifstream input(annotation_file_name.c_str());
    if(!input.is_open()) return;

    Chain* chain = (protein->chains[0]);
    int residues = chain->getResidues()[chain->getResidues().size()-1]->getId()+1;//Index of last residue + 1
    protein->residueAnnotations = new int[residues];
    for(int i=0;i<residues;i++)
        protein->residueAnnotations[i] = 0;

    int residue, annotation;
    while (input.good()) {
        input>>residue;
        input>>annotation;
        if(residue<0){
            cerr<<"IO::readAnnotations: Error - Residues can not have negative IDs"<<endl;
            exit(-1);
        }
        if(residue>=residues){
            cerr<<"IO::readAnnotations: Error - Trying to annotate residue "<<residue<<" but largest res is "<<(residues-1)<<endl;
            exit(-1);
        }

        assert(residue<residues);
        protein->residueAnnotations[residue] = annotation;
        //log("debugRas")<<"resAnn["<<residue<<"] = "<<annotation<<endl;
    }

}

void IO::readHbonds_dssr(Molecule * protein, string dssrFile){
	cerr<<"IO::readHbonds_dssr(..) not yet implemented. Sorry"<<endl;
	exit(-1);
}

void IO::readHbonds_rnaview(Molecule * protein, string file, bool fillAnnotations){
  int residues=0;
  if(fillAnnotations){
    for(auto const& chain: protein->chains){
      int lastRes = chain->getResidues()[chain->getResidues().size()-1]->getId()+1;//Index of last residue + 1
      if(lastRes>residues) residues = lastRes;
    }
    protein->residueAnnotations = new int[residues];
    for(int i=0;i<residues;i++)
      protein->residueAnnotations[i] = 0;
  }


  ifstream input(file.c_str());
  if(!input.is_open()) {cerr<<"IO::readHbonds_rnaview(..) Error: Couldnt open "<<file<<" for reading"<<endl; exit(-1);}

  string line; bool readingBasePairs = false;
  while( getline(input, line) ){
    if(line.at(line.length()-1)=='\r') line = line.substr(0,line.length()-1);
    if(line=="BEGIN_base-pair"){ readingBasePairs = true; continue; }
    if(line=="END_base-pair") break;
    if(!readingBasePairs) continue;

    string chain1 = line.substr(11,1);
    string chain2 = line.substr(30,1);
    cerr<<"IO::readHbonds_rnaview - TODO: check that chains are read correctly ("<<chain1<<", "<<chain2<<")"<<endl;
    int res1 = atoi(Util::trim(line.substr(13,6)).c_str());
    int res2 = atoi(Util::trim(line.substr(23,6)).c_str());

    if(fillAnnotations){
      if(res1>0 && res1<residues)	protein->residueAnnotations[res1] = 1;
      if(res2>0 && res2<residues)	protein->residueAnnotations[res2] = 1;
    }

    if(line.length()<50) continue;
    string type = Util::trim(line.substr(49));
    if(type.at(type.length()-1)=='\r') type = type.substr(0,type.length()-1);

    string donorName, acceptorName;
    if(type=="XIX"){ //Watson-Crick G-C
      if(line.at(20)=='C'){ int tmp=res1; res1=res2; res2=tmp; }//Swap. res1 (donor) is now G and res2 (acceptor) is C
      donorName = "N1";
      acceptorName = "N3";
    }else if(type=="XX"){ //Watson-Crick A-U
      if(line.at(20)=='A'){ int tmp=res1; res1=res2; res2=tmp; }//Swap. res1 (donor) is now U and res2 (acceptor) is A
      donorName = "N3";
      acceptorName = "N1";
    }else if(type=="XXI"){ //Protonated A-C (selfmade type)
      if(line.at(20)=='C'){ int tmp=res1; res1=res2; res2=tmp; }//Swap. res1 (donor) is now A and res2 (acceptor) is C
      donorName = "N1";
      acceptorName = "O2"; //(This might primarily be valid for 1LDZ)
    }else continue;
    Atom* donor = protein->getAtom(chain1, res1, donorName);
    if(donor==NULL) { cerr<<"IO::readHbonds_rnaview(..) Error: Couldnt find "<<res1<<"_"<<donorName<<endl;exit(-1);}
    Atom* acceptor = protein->getAtom(chain2, res2, acceptorName);
    if(acceptor==NULL) { cerr<<"IO::readHbonds_rnaview(..) Error: Couldnt find "<<res2<<"_"<<acceptorName<<endl;exit(-1);}
    Atom* AA = acceptor->getFirstCovNeighbor();
    Atom* hatom = NULL;
    for(vector<Atom*>::iterator ait = donor->Cov_neighbor_list.begin(); ait != donor->Cov_neighbor_list.end(); ait++){
      Atom* a = *ait;
      if(a->Element==atomH){ hatom = a; break; }
    }
    if(hatom==NULL){//TODO: Manually create H-atom
      //int newId = *(m_protein->atoms.end()-1)
      //Atom* atom = new Atom("H",,pos,occ,B);
      //m_protein->addAtom(chain_name,res_name,res_id,atom);

      cerr<<"IO::readHbonds_rnaview(..). Error: Couldnt locate hydrogen in residue "<<res1<<endl;
      cerr<<"IO::readHbonds_rnaview(..). Error: "<<line<<endl;
      exit(-1);
    }
    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    protein->addHbond(new_hb);
    log("so")<<"RNAView file: Created hydrogen bond between "<<hatom<<
        "("<<hatom->getId()<<") and "<<acceptor<<"("<<acceptor->getId()<<")"<<endl;
  }

  if(fillAnnotations){
    for(int i=0;protein->residueAnnotations[i]==0;i++)
      protein->residueAnnotations[i]=1;
    for(int i=residues-1;protein->residueAnnotations[i]==0;i--)
      protein->residueAnnotations[i]=1;
  }
  //enableLogger("debugRas");
  log("so")<<"ResidueAnnotations: ";
  for(int i=0;i<residues;i++){
    log("so")<<protein->residueAnnotations[i];
  }
  log("so")<<endl;

}

void IO::readHbonds_first(Molecule * protein, string file){
  ifstream input(file.c_str());
  if(!input.is_open()) {cerr<<"IO::readHbonds_first(..) Error: Couldnt open "<<file<<" for reading"<<endl; exit(-1);}

  string line;
  while( getline(input, line) ){
    if(line.at(line.length()-1)=='\r') line = line.substr(0,line.length()-1);

    int id1 = atoi(Util::trim(line.substr(0,8)).c_str());
    int id2 = atoi(Util::trim(line.substr(8,8)).c_str());

    Atom* hatom = protein->getAtom(id1);
    Atom* acceptor = protein->getAtom(id2);
    Atom* AA = acceptor->getFirstCovNeighbor();
    Atom* donor = hatom->getFirstCovNeighbor();

    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    protein->addHbond(new_hb);
    log("verbose")<<"FIRST file: Created hydrogen bond between "<<hatom<<" and "<<acceptor<<endl;
  }
}


void IO::readHbonds_vadar(Molecule * protein, string file){
  ifstream input(file.c_str());
  if(!input.is_open()) {cerr<<"IO::readHbonds_vadar(..) Error: Couldnt open "<<file<<" for reading"<<endl; exit(-1);}

  string line;
  string res1, name1, res2, name2;
  double dist;
  bool reading = true;
  while( reading ){
    //if(line.at(line.length()-1)=='\r') line = line.substr(0,line.length()-1);
    //reading = input>>res1>>name1>>res2>>name2>>dist;
    input>>res1>>name1>>res2>>name2>>dist;
    reading = input.good();
    //cout<<"Res1: "<<res1<<" Name1: "<<name1<<" Res2: "<<res2<<" Name2: "<<name2<<" dist: "<<dist<<endl ;

    string chain1 = res1.substr(res1.length()-1,1);
    string chain2 = res2.substr(res2.length()-1,1);
    int resId1 = atoi(res1.substr(0,res1.length()-1).c_str());
    int resId2 = atoi(res2.substr(0,res2.length()-1).c_str());
    //TODO: Assumes the chains are the same

    Atom* donor = protein->getAtom(chain1, resId1, name1);
    if(donor==NULL) { cerr<<"Cannot find res "<<resId1<<"/"<<name1<<endl; exit(-1); }
    Atom* acceptor = protein->getAtom(chain1, resId2, name2);
    if(acceptor==NULL) { cerr<<"Cannot find res "<<resId2<<"/"<<name2<<endl; exit(-1); }
    cout<<"Hbond "<<donor->getId()<<" "<<acceptor->getId()<<endl;
    Atom* hatom = donor->getFirstCovNeighbor();
    Atom* AA = acceptor->getFirstCovNeighbor();
    //Atom* hatom = donor->getFirstCovHydrogenNeighbor();
    //Atom* AA = acceptor->getFirstCovNeighbor();
    //if(hatom==NULL){
    //	hatom = donor->getFirstCovHydrogenNeighbor();
    //	h
    //}
    //Atom* AA = hatom->getFirstCovNeighbor();

    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    protein->addHbond(new_hb);
    log("verbose")<<"vadar file: Created hydrogen bond between "<<hatom<<" and "<<acceptor<<endl;
  }

}

void IO::writePyMolScript(Molecule * protein, string pdb_file, string output_file_name) {

	ofstream pymol_script( output_file_name.c_str() );

	string render_style;

	  ostringstream flexible;
	  flexible << "create Flexible, ";

	  pymol_script << "# Rigid cluster coloring script for PyMol" << endl << "#" << endl;
	  pymol_script << "# Original by Dan Farrell, Brandon Hespenheide." << endl;
	  pymol_script << "# Arizona State University" << endl;
	  pymol_script << "# Adapted by Dominik Budday." << endl;
	  pymol_script << "# Friedrich-Alexander University of Erlangen-Nuremberg" << endl;
	  pymol_script << "############################################################" << endl << endl;

	  pymol_script << "from pymol import cmd" << endl;
	  pymol_script << "from pymol.cgo import *" << endl << endl;

	  // Some final global attributes to set.
	  //////////////////////////////////////////////////////////////////////
	  pymol_script << "bg_color white" << endl;

	  // Color the whole structure black with thin lines. This represents the
	  // underlying flexible network.
	  //////////////////////////////////////////////////////////////////////
      pymol_script << "load " << pdb_file << endl << endl;
	  pymol_script << "set line_width = 1" << endl;
	  pymol_script << "color black" << endl << endl;

	  // For each rigid cluster, calculate a unique color, and have pymol
	  // assign that color to those atoms.
	  // Also, we save each cluster's color in a map for use again later
	  map< unsigned int, string > mapClusterIDtoColor;
	  map<unsigned int,Rigidbody*>::iterator rbit = protein->m_conf->m_biggerRBMap.begin();
	  vector< pair< int, unsigned int> >::iterator sit = protein->m_conf->m_sortedRBs.begin();
	  Color::next_rcd_color_as_name(true);

//	  while(rbit != m_protein->m_conf->m_biggerRBMap.end() ){
//		  if( (rbit->second)->size() >= MIN_CLUSTER_SIZE ){
//			  string color = Color::next_rcd_color_as_name();
//			  mapClusterIDtoColor[rbit->first] = color;
//			  pymol_script << "color " << color << ", ( b > " << float(rbit->first-0.01)
//			  << " and b < " << float(rbit->first+0.01) << ")" << endl;
//		  }
//		  rbit++;
//	  }

	  while(sit != protein->m_conf->m_sortedRBs.end() ){
	  		  if( sit->first >= MIN_CLUSTER_SIZE ){
//	  			  cout<<"Current print: "<<sit->second<<endl;
	  			  string color = Color::next_rcd_color_as_name();
	  			  mapClusterIDtoColor[sit->second] = color;
	  			  pymol_script << "color " << color << ", ( b > " << float(sit->second-0.01)
	  			  << " and b < " << float(sit->second+0.01) << ")" << endl;
	  		  }
	  		  sit++;
	  }

	  // Python commands to draw hbonds
	    //////////////////////////////////////////////////////////////////////
	      pymol_script << "# Draw hbonds as distance objects" << endl;
	      pymol_script << "set dash_gap, 0.1" << endl;

	      int site_1, site_2;
	      vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator eit;

	      for (eit=protein->m_spanning_tree->CycleAnchorEdges.begin(); eit != protein->m_spanning_tree->CycleAnchorEdges.end(); eit++) {

	          site_1 = eit->first->getBond()->Atom1->getId();
	          site_2 = eit->first->getBond()->Atom2->getId();
	          pymol_script << "distance hbonds = id " << site_1 << " , id " << site_2 << endl;

	      }
	      pymol_script << "color red, hbonds" << endl;
	      pymol_script << "hide labels, hbonds" << endl;

	  // Create the rigid cluster objects for pymol, and color them. Only those
	  // clusters larger than the min_output_cluster_size will have objects
	  // created for them. The remaining clusters will be binned into bulk
	  // objects.
	  //////////////////////////////////////////////////////////////////////

	  render_style = "lines";

	  // for each rigid cluster, create a new pymol object of the cluster.
	  unsigned int total_RC_objects = 0;
	  rbit = protein->m_conf->m_biggerRBMap.begin();
	  sit = protein->m_conf->m_sortedRBs.begin();
	  while( sit != protein->m_conf->m_sortedRBs.end() ){
		  if( (sit->first) >= MIN_CLUSTER_SIZE ){
		  			  pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << sit->second<< " has "
		  			  << sit->first << " atoms." << endl;
		  			  pymol_script << "create RC" << sit->second << ", ( b > " << float(sit->second-0.01)
		  			  << " and b < " << float(sit->second+0.01) << ")" << endl;
		  			  pymol_script << "show " << render_style << ", RC" << sit->second << endl;
		  			  pymol_script << "set line_width = 3, " << "RC" << sit->second << endl;
		  			  pymol_script << "color " << mapClusterIDtoColor[sit->second] << ", RC" << sit->second << endl << endl;
		  }
		   sit++;
	  }

	  // Some final global attributes to set.
	  //////////////////////////////////////////////////////////////////////
	  pymol_script << "bg_color white" << endl;
	  pymol_script << "clip slab, 200" << endl;
	  pymol_script << "center all" << endl;
	  pymol_script << "zoom" << endl <<endl;

	  pymol_script << "hide everything" << endl;
	  pymol_script << "show lines" << endl;
	  pymol_script << "show cartoon" << endl;
	  pymol_script << "show dashes, hbonds" << endl;
	  pymol_script << "set cartoon_fancy_helices, on" << endl;
	  pymol_script << "set ray_opaque_background, on" << endl;


	  pymol_script.close();
}


void IO::writeRBs(Molecule * protein, string output_file_name){

	ofstream output( output_file_name.c_str() );
	if(!output.is_open()) {
			cerr<<"Cannot write to "<<output_file_name<<". You might need to create output directory first"<<endl;
			exit(-1);
	}
	if(protein->m_conf!=NULL){
		Configuration* c = protein->m_conf;
		vector< pair<int, unsigned int> >::iterator it=c->m_sortedRBs.begin();

		while( it != c->m_sortedRBs.end() ){

			vector<Atom*>::iterator ait = c->m_biggerRBMap[it->second]->Atoms.begin();
			while( ait != c->m_biggerRBMap[it->second]->Atoms.end() ){
				output << (*ait)->getId() << endl;
				ait++;
			}
			output << endl;
			output << "NaN"<<endl;
			output << endl;
			it++;
		}

	}
	output.close();

}


void IO::writeStats(Molecule * protein, string output_file_name){

	ofstream output( output_file_name.c_str() );
	if(!output.is_open()) {
			cerr<<"Cannot write to "<<output_file_name<<". You might need to create output directory first"<<endl;
			exit(-1);
	}
	if(protein->m_conf!=NULL){
		int diff= protein->m_spanning_tree->num_DOFs - protein->m_spanning_tree->Cycle_DOF_num;
		int sum = protein->m_conf->getNullspace()->NullspaceSize() + diff;
		//int sum = m_protein->m_conf->CycleNullSpace->getNullspace()Size + diff;

		output <<"************* Statistics on the molecular structure *************"<<endl;
		output << "Number of chains: " << protein->chains.size() << endl;
		output << "Number of atoms: " << protein->atoms.size() << endl;
		output <<"Number of covalent bonds: " << protein->Cov_bonds.size()<< endl;
		output <<"Number of hydrogen bonds: " << protein->m_spanning_tree->CycleAnchorEdges.size()<< endl;
		output << "Number of dihedrals in spanning tree: " << protein->m_spanning_tree->num_DOFs << endl;
		output <<"Number of free dihedrals: " << diff << endl;
		output <<"Number of cycle dihedrals: " << protein->m_spanning_tree->Cycle_DOF_num << endl<<endl;
		output <<"************* Statistics on rigidity analysis *************"<<endl;
		output <<"Number of internal m_dofs (nullspace dimension): " << protein->m_conf->getNullspace()->NullspaceSize() << endl;
		//output <<"Number of internal m_dofs (nullspace dimension): " << m_protein->m_conf->CycleNullSpace->getNullspace()Size <<endl;
		output << "Overall number of m_dofs (free + coordinated): " << sum << endl;
		output <<"Number of rigidified covalent bonds: " << protein->m_conf->getNullspace()->NumRigidDihedrals() << endl;
		//output <<"Number of rigidified covalent bonds: " << m_protein->m_conf->CycleNullSpace->m_numRigid<<endl;
		output <<"Number of rigidified hydrogen bonds: " << protein->m_conf->getNullspace()->NumRigidHBonds() << endl;
		//output <<"Number of rigidified hydrogen bonds: " << m_protein->m_conf->CycleNullSpace->m_numRigidHBonds<<endl;
		output << "Number of rigid clusters: " << protein->m_conf->m_numClusters << endl;
		output << "Size of biggest cluster: " << protein->m_conf->m_maxSize << endl;
		output << "Index of biggest cluster: " << protein->m_conf->m_maxIndex<< endl;
	}
	output.close();

}

void IO::writeTrajectory (Molecule*molecule, string output_file_name, string output_mdl, Molecule* target) {

	ofstream output(output_file_name.c_str());
	if(!output.is_open()) {
		//cerr<<"Cannot write to "<<output_file_name<<". You might need to create output directory first"<<endl;
		cerr<<"Cannot write to "<<output_file_name<<endl;
		exit(-1);
	}
	int firstPathLength;
	Configuration* cFinalFwd = molecule->m_conf;

	if(molecule->m_conf!=NULL){

		Configuration* c = molecule->m_conf;
		int currId = c->m_id;
		int treeDepth = c->m_treeDepth;
		firstPathLength=treeDepth;

		output << "REMARK\tFound connected path of "<<treeDepth<<" samples."<<endl;
		output << "REMARK\tFinal sample of first path is "<<currId<<endl;
		output << "REMARK\tDistance to target is "<<setprecision(6)<<c->m_distanceToTarget<<endl;
		output << "REMARK\tPath is "<<currId;

		//Make sure to follow the tree from initial to final
		double *path = new double[treeDepth+1];
		path[treeDepth]=currId;
		for(int i=treeDepth-1; i>=0; i--){
			c = c->getParent();
			currId = c->m_id;
			output<<", "<<currId;
			path[i]=currId;
		}
		output<<endl;

		//Align configuration!
//		Rigidbody* biggestRb = c->m_biggerRBMap.
//		list<unsigned int> *atom_ids = c->m_sortedRBs;
		//Set initial configuration
//		molecule->SetConfiguration(c);
		c->updateProtein();
//		molecule->alignToReference()
		c = molecule->m_conf;
		currId = c->m_id;
		Configuration* iniConf = c;
		output << "REMARK\tInitial Distance to target was "<<setprecision(6)<<c->m_distanceToTarget<<endl;

		vector<Edge*>::iterator eit;

		for(int i=0; i <= treeDepth; i++){

			output << "MODEL"<<currId<<endl;

//			std::size_t pos = output_file_name.find("path");
//			string input_file_name=output_file_name.substr(0,pos) + "new_" +
//					static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".pdb";
//			ifstream inputFile(input_file_name.c_str());
//
//			if (!inputFile.good()) {
//				cerr << "Error: cannot open file " << input_file_name << endl;
//				exit(-1);
//			}
//			string line, temp;
//
//			int test=0;
//			while (!inputFile.eof()) {
//				getline(inputFile,line);
//				if (inputFile.eof()) break;
//				if (line.substr(0,4)!="ATOM") // skip if it is not an ATOM line
//					continue;
//				output << line << endl;
//				test++;
//			}

//			map<unsigned int, unsigned int> resiColorMap;
//			for( eit = molecule->m_spanning_tree->Edges.begin(); eit != molecule->m_spanning_tree->Edges.end(); eit++){
//				Edge* e = (*eit);
//				int dofId = e->DOF_id;
//				CTKResidue* res = e->Bond->Atom1->Parent_residue;
//				double val = Abs(c->m_f[dofId]);
//				int num;
//				if( val > 0.01)
//					num=1;
//				else if (val > 0.00000001)
//					num=2;
//				else
//					num=3;
//				if( resiColorMap.find(res->Id) != resiColorMap.end()){
//					if( num < resiColorMap.find(res->Id)->second)
//						resiColorMap.insert( make_pair(res->Id, num) );
//				}
//				else{
//					resiColorMap.insert( make_pair(res->Id, num) );
//				}
//			}

			for (auto const& atom: molecule->atoms) {
				Residue* res = atom->getResidue();
				char buffer[100];
				sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6d          %2s  ",
								atom->getId(),atom->getName().c_str(),
								res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
						//atom->Position.x,atom->Position.y,atom->Position.z,atom->getType().c_str());
								atom->m_Position.x,atom->m_Position.y,atom->m_Position.z,atom->getBiggerRigidbody()->id(),atom->getType().c_str() );
//					atom->Position.x,atom->Position.y,atom->Position.z,resiColorMap.find(res->getId())->second,atom->getType().c_str() );
				string line(buffer);
				output << line << endl;
			}

			output << "ENDMDL"<<endl;
//			string out_q = output_file_name.substr(0,output_file_name.length()-8) + "q_"+ static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".txt";
//			writeQ(molecule, iniConf, out_q);
//			string out_pdb = output_file_name.substr(0,output_file_name.length()-8) + "new_" + static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".pdb";
//			writePdb(molecule, out_pdb);

			//Next configuration only if we haven't reached the last one yet!
			if( i < treeDepth){

				ConfigurationList::iterator cit;
				for(auto const& cit: c->getChildren()){
					if(cit->m_id != path[i+1])
						continue;
					else{
//						molecule->SetConfiguration(*cit);
						cit->updateProtein();
						break;
					}
				}
				c = molecule->m_conf;
				currId = c->m_id;
			}
		}
	}
	if( target != NULL && target->m_conf != NULL){//also add the additional path of the reversePlanner

		Configuration* c = target->m_conf;
		int currId = c->m_id;
		int treeDepth = c->m_treeDepth;

		double alignVal = metrics::RMSD::align(molecule,target);
		double distance = metrics::RMSD::distance_noOptimization(cFinalFwd,c);

		output << "REMARK\tAdding reversePlanning path of "<<treeDepth<<" samples."<<endl;
		output << "REMARK\tOverall path length is "<<firstPathLength+treeDepth<<" samples."<<endl;
		output << "REMARK\tStarting from reverse sample "<<currId<<endl;
		output << "REMARK\tRMSD at transition is "<<setprecision(3)<<distance<<endl;
		output << "REMARK\tReverse path follows "<<currId;

		while( c->getParent() != NULL){
			c=c->getParent();
			currId = c->m_id;
			output << ", "<<currId;
		}
		output <<endl;
		c = target->m_conf;
		currId = c->m_id;
		Configuration* iniConf = c;

		vector<Edge*>::iterator eit;

		for(int i=0; i <= treeDepth; i++){

			output << "REMARK\tFile corresponds to sample "<<currId<<endl;
			output << "MODEL"<<++firstPathLength<<endl;

//			std::size_t pos = output_file_name.find("path");
//			string input_file_name=output_file_name.substr(0,pos) + "new_" +
//					static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".pdb";
//			ifstream inputFile(input_file_name.c_str());
//
//			if (!inputFile.good()) {
//				cerr << "Error: cannot open file " << input_file_name << endl;
//				exit(-1);
//			}
//			string line, temp;
//
//			int test=0;
//			while (!inputFile.eof()) {
//				getline(inputFile,line);
//				if (inputFile.eof()) break;
//				if (line.substr(0,4)!="ATOM") // skip if it is not an ATOM line
//					continue;
//				output << line << endl;
//				test++;
//			}
//			map<unsigned int, unsigned int> resiColorMap;
//			for( eit = target->m_spanning_tree->Edges.begin(); eit != target->m_spanning_tree->Edges.end(); eit++){
//				Edge* e = (*eit);
//				int dofId = e->DOF_id;
//				CTKResidue* res = e->Bond->Atom1->Parent_residue;
//				double val = Abs(c->m_f[dofId]);
//				int num;
//				if( val > 0.01)
//					num=1;
//				else if (val > 0.00000001)
//					num=2;
//				else
//					num=3;
//				if( resiColorMap.find(res->Id) != resiColorMap.end()){
//					if( num < resiColorMap.find(res->Id)->second)
//						resiColorMap.insert( make_pair(res->Id, num) );
//				}
//				else{
//					resiColorMap.insert( make_pair(res->Id, num) );
//				}
//			}

			for (auto const& atom: target->atoms) {
				Residue* res = atom->getResidue();
				char buffer[100];
				//sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  ",
				sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6d          %2s  ",
								atom->getId(),atom->getName().c_str(),
								res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
						//atom->Position.x,atom->Position.y,atom->Position.z,atom->getType().c_str());
								atom->m_Position.x,atom->m_Position.y,atom->m_Position.z,atom->getBiggerRigidbody()->id(),atom->getType().c_str() );
//					atom->Position.x,atom->Position.y,atom->Position.z,resiColorMap.find(res->getId())->second,atom->getType().c_str() );
				string line(buffer);
				output << line << endl;
			}

			output << "ENDMDL"<<endl;
//			string out_q = output_file_name.substr(0,output_file_name.length()-8) + "rev_q_"+ static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".txt";
//			writeQ(target, iniConf, out_q);
//			string out_pdb = output_file_name.substr(0,output_file_name.length()-8) + "rev_new_" + static_cast<ostringstream*>( &(ostringstream() << currId) )->str() + ".pdb";
//			writePdb(target, out_pdb);

			//Next configuration only if we haven't reached the last one yet!
			if( c->getParent() != NULL){
				c=c->getParent();
//				target->SetConfiguration(c);
				c->updateProtein();
				currId = c->m_id;
			}
		}
	}

	output.close();

	ofstream pymol_script( output_mdl.c_str() );

	string render_style;
//
//	  ostringstream flexible;
//	  flexible << "create Flexible, ";

	pymol_script << "# Colored Movie for PyMol" << endl << "#" << endl;
	pymol_script << "# Inspired by Dan Farrell, Brandon Hespenheide." << endl;
	pymol_script << "# Author: Dominik Budday." << endl;
	pymol_script << "# Friedrich-Alexander University of Erlangen-Nuremberg" << endl;
	pymol_script << "############################################################" << endl << endl;

	pymol_script << "from pymol import cmd" << endl;
	pymol_script << "from pymol.cgo import *" << endl << endl;

	// Some final global attributes to set.
	//////////////////////////////////////////////////////////////////////
	pymol_script << "bg_color white" << endl;

	// Color the whole structure black with thin lines. This represents the
	// underlying flexible network.
	//////////////////////////////////////////////////////////////////////
	pymol_script << "load " << output_file_name << endl << endl;
	pymol_script << "set line_width = 1" << endl;
//	  pymol_script << "color black" << endl << endl;

	//Determine object name
	size_t start_idx = output_file_name.find_last_of("/")+1;
	size_t end_idx = output_file_name.find_last_of(".");
	string obj_name = output_file_name.substr(start_idx, end_idx-start_idx);
	pymol_script << "intra_fit " << obj_name << endl << endl;


	// For each rigid cluster, calculate a unique color, and have pymol
	// assign that color to those atoms.
	// Also, we save each cluster's color in a map for use again later

	map< unsigned int, string > mapClusterIDtoColor;
	map<unsigned int,Rigidbody*>::iterator rbit = molecule->m_conf->m_biggerRBMap.begin();
	vector< pair< int, unsigned int> >::iterator sit = molecule->m_conf->m_sortedRBs.begin();
	Color::next_rcd_color_as_name(true);

	while(sit != molecule->m_conf->m_sortedRBs.end() ){
		if( sit->first >= MIN_CLUSTER_SIZE ){
			string color = Color::next_rcd_color_as_name();
			mapClusterIDtoColor[sit->second] = color;
			pymol_script << "color " << color << ", ( b > " << float(sit->second-0.01)
			<< " and b < " << float(sit->second+0.01) << ")" << endl;
		}
		sit++;
	}

	// Python commands to draw hbonds
	//////////////////////////////////////////////////////////////////////
	pymol_script << "# Draw hbonds as distance objects" << endl;
	pymol_script << "set dash_gap, 0.1" << endl;

	int site_1, site_2;
	vector< pair<Edge*,RigidbodyGraphVertex*> >::iterator eit;

	for (auto  const& eit: molecule->m_spanning_tree->CycleAnchorEdges){

		site_1 = eit.first->getBond()->Atom1->getId();
		site_2 = eit.first->getBond()->Atom2->getId();
		pymol_script << "distance hbonds = id " << site_1 << " , id " << site_2 << endl;

	}
	pymol_script << "color pink, hbonds" << endl;
	pymol_script << "hide labels, hbonds" << endl;

	// Create the rigid cluster objects for pymol, and color them. Only those
	// clusters larger than the min_output_cluster_size will have objects
	// created for them. The remaining clusters will be binned into bulk
	// objects.
	//////////////////////////////////////////////////////////////////////

	render_style = "lines";
//	   for each rigid cluster, create a new pymol object of the cluster.
	unsigned int total_RC_objects = 0;
	rbit = molecule->m_conf->m_biggerRBMap.begin();
	sit = molecule->m_conf->m_sortedRBs.begin();
	while( sit != molecule->m_conf->m_sortedRBs.end() ){
		if( (sit->first) >= MIN_CLUSTER_SIZE ){
			pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << sit->second<< " has "
			<< sit->first << " atoms." << endl;
			pymol_script << "create RC" << sit->second << ", ( b > " << float(sit->second-0.01)
			<< " and b < " << float(sit->second+0.01) << ")" << endl;
			pymol_script << "show " << render_style << ", RC" << sit->second << endl;
			pymol_script << "set line_width = 3, " << "RC" << sit->second << endl;
			pymol_script << "color " << mapClusterIDtoColor[sit->second] << ", RC" << sit->second << endl << endl;
		}
		sit++;
	}

//	  pymol_script << "create Moving, ( b > " << float(0.9)
//		<< " and b < " << float(1.1) << ")" << endl;
//	  pymol_script << "color red, Moving"<< endl;
//
//	  pymol_script << "create Adaptive, ( b > " << float(1.9)
//		<< " and b < " << float(2.1) << ")" << endl;
//	  pymol_script << "color yellow, Adaptive"<< endl;
//
//	  pymol_script << "create Fix, ( b > " << float(2.9)
//		<< " and b < " << float(3.1) << ")" << endl;
//	  pymol_script << "color blue, Fix"<< endl;

	// Some final global attributes to set.
	//////////////////////////////////////////////////////////////////////
	pymol_script << "bg_color white" << endl;
	pymol_script << "clip slab, 200" << endl;
	pymol_script << "center all" << endl;
	pymol_script << "zoom" << endl <<endl;

	pymol_script << "hide everything" << endl;
	pymol_script << "show lines" << endl;
	pymol_script << "set line_width = 3 " << endl;
	pymol_script << "show cartoon" << endl;
	pymol_script << "show dashes, hbonds" << endl;
	pymol_script << "set cartoon_transparency, 0.7" << endl;

	pymol_script.close();
}