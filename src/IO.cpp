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
#include <math/NullspaceSVD.h>


#include "CTKTimer.h"
#include "IO.h"
#include "core/Atom.h"
#include "core/Chain.h"
#include "core/Bond.h"
#include "core/HBond.h"
#include "DisjointSets.h"
#include "Logger.h"
#include "Selection.h"
#include "Color.h"
#include "HbondIdentifier.h"


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

Molecule* IO::readPdb (
    const string& pdb_file,
    Selection movingResidues,
    const vector<string>& extraCovBonds,
    const vector<int>& roots,
    const string& hbondMethod,
    const string& hbondFile,
    const Molecule* reference
)
{
  Molecule* molecule = new Molecule();

  //Read molecule name
  int nameSplit=pdb_file.find_last_of("/\\");
  string moleculeName=pdb_file.substr(nameSplit + 1);
  int pos=moleculeName.rfind(".pdb");
  if( pos != string::npos ) moleculeName=moleculeName.substr(0, pos);
  molecule->setName(moleculeName);

  //Open file
  ifstream pdb(pdb_file.c_str());
  if( !pdb.good()) {
    cerr << "Error: cannot open file " << pdb_file << endl;
    exit(1);
  }
  string line, temp;
  vector< vector<int > > conectRecords;
  vector< pair<int,int> > torsionConstraints;

  int atomCount=0;

  //Read lines of PDB-file.
  while( !pdb.eof()) {
    getline(pdb, line);
    if( pdb.eof()) break;

    //Parse CONECT records for extra covalent bonds (typically ligands)
    if( line.substr(0, 6) == "CONECT" ){
      vector<int> conectEntry;
      string remainder = line.substr(6);
      int numRecs = int(remainder.length()/5);
      for(int rec=0; rec<numRecs;rec++){
        int entry = atoi(remainder.substr(rec*5,5).c_str());
        conectEntry.push_back(entry);
      }
      conectRecords.push_back(conectEntry);
      continue;
    }

    //Parse REMARK records
    if( line.substr(0,29)=="REMARK 555 RevoluteConstraint"){
      std::istringstream iss(line.substr(29));
      int id1, id2;
      //TODO: Fail gracefully
      iss>>id1>>id2;
      torsionConstraints.push_back( make_pair(id1,id2) );

      continue;
    }

    //Parse ATOM/HETATM records and add them to molecule
    if( line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM" ) {
      ++atomCount;
      int offset = 0;
      if (atomCount >99999) //quick and dirty adaption for molecules with > 100000 atoms
        offset = 1;
      // chain info
      string chain_name = line.substr(21+offset, 1); // line[22]
      // residue info
      int res_id = atoi(line.substr(22+offset, 4).c_str()); // line[23:26]
      string res_name = Util::trim(line.substr(17+offset, 3)); // line[18:20]
      // atom info
      int atom_id = atoi(line.substr(6+offset, 5).c_str()); // line[7:11]s
      string atom_name = Util::trim(line.substr(12+offset, 5)); // line[13:17]
      if (atom_name == "OP3") continue;
      if (atom_name.at(0) >= 49 && atom_name.at(0) <= 57) { // if the first char is 1-9
        string temp_name(atom_name.substr(1, 3));
        temp_name += atom_name.substr(0, 1);
        atom_name = temp_name;
      }

      double x = atof(line.substr(30+offset, 8).c_str()); // line[31:38]
      double y = atof(line.substr(38+offset, 8).c_str()); // line[39:46]
      double z = atof(line.substr(46+offset, 8).c_str()); // line[47:54]

      Coordinate pos(x, y, z);
      Atom* atom = molecule->addAtom(chain_name, res_name, res_id, atom_name, atom_id, pos);

      float occupancy = 0.0;
      if(line.length()>=59) occupancy = atof(line.substr(56+offset,4).c_str());
      atom->setOccupancy(occupancy);

      float bfactor = 0.0;
      if(line.length()>=65) bfactor = atof(line.substr(60+offset,6).c_str());
      atom->setBFactor(bfactor);

      continue;
    }
  }
  pdb.close();

  //Warn if no atoms read
  if(molecule->getAtoms().size()==0){
    cerr<<"IO::readPdb - No atoms read from file "<<pdb_file<<endl;
    exit(-1);
  }

  //Warn if no hydrogens found
  bool foundHydro = false;
  for(auto const& atom: molecule->getAtoms()){
    if(atom->m_element==AtomType::atomH){
      foundHydro = true;
      break;
    }
  }
  if(!foundHydro) cerr<<"IO::readPdb - Warning: No hydrogens found in file. Consider running `reduce`"<<endl;


  //Create covalent bonds from residue profiles
  ResidueProfile residue_profile=readResidueProfile();
  for( auto const &cur_chain: molecule->chains ) {
    for( auto const &cur_res: cur_chain->getResidues()) {
      string res_name=cur_res->getName();
      map<string, vector<CovBond> >::iterator profile_it=residue_profile.find(res_name);

      if( profile_it == residue_profile.end()) {
        cerr << "IO::readPdb - warning: Unknown residue " << res_name << ". Atoms will have fixed positions." << endl;
        continue;
      }

      for( auto const& profile_bond: profile_it->second ){
        string atom_name1 = profile_bond.first;
        string atom_name2 = profile_bond.second;
        Residue *res1 = cur_res;
        Residue *res2 = cur_res;
        if( atom_name1[0] == '-' ) {
          atom_name1=Util::trim(atom_name1, '-');
          res1=cur_res->getPrevResidue();
        }
        else if( atom_name1.at(0) == '+' ) {
          atom_name1=Util::trim(atom_name1, '+');
          res1=cur_res->getNextResidue();
        }
        if( atom_name2.at(0) == '-' ) {
          atom_name2=Util::trim(atom_name2, '-');
          res2=cur_res->getPrevResidue();
        }
        else if( atom_name2.at(0) == '+' ) {
          atom_name2=Util::trim(atom_name2, '+');
          res2=cur_res->getNextResidue();
        }
        if( res1 == nullptr || res2 == nullptr ) {
          continue;
        }

//        makeCovBond(res1, res2, atom_name1, atom_name2);
//        molecule->addCovBond(res1,res2,atom_name1,atom_name2);
        Atom* atom1 = res1->getAtom(atom_name1);
        Atom* atom2 = res2->getAtom(atom_name2);
        if(atom1 && atom2)
          molecule->addCovBond(atom1, atom2);
      } // finish looping over bonds
    } // finish looping over residues
  } // finish looping over chains

  //Create covalent bonds from extraCovBonds variable
  if( reference == nullptr ) {
    for( unsigned int i=0; i<extraCovBonds.size(); i++ ) {
      vector<string> tokens=Util::split(extraCovBonds[i], '-');
      int atomId1=atoi(tokens[0].c_str());
      int atomId2=atoi(tokens[1].c_str());
      Atom *a1=molecule->getAtom(atomId1);
      Atom *a2=molecule->getAtom(atomId2);
      if( a1 == nullptr ) {
        cerr << "Cannot find atom with id " << atomId1 << endl;
        exit(-1);
      }
      if( a2 == nullptr ) {
        cerr << "Cannot find atom with id " << atomId2 << endl;
        exit(-1);
      }
//      makeCovBond(a1->getResidue(), a2->getResidue(), a1->getName(), a2->getName());
//      molecule->addCovBond(a1->getResidue(), a2->getResidue(),a1->getName(), a2->getName());
      molecule->addCovBond(a1,a2);
      log("dominik") << "Creating bond between " << a1 << " and " << a2 << " in protein " << molecule->getName() << endl;
    }
  }
  else {
    for( unsigned int i=0; i<extraCovBonds.size(); i++ ) {
      vector<string> tokens=Util::split(extraCovBonds[i], '-');
      int atomId1=atoi(tokens[0].c_str());
      int atomId2=atoi(tokens[1].c_str());
      Atom *a1=reference->getAtom(atomId1);//id's from reference
      Atom *a2=reference->getAtom(atomId2);
      //residue Ids and names
      int resId1=a1->getResidue()->getId();
      string name1=a1->getName();
      string chainName1=a1->getResidue()->getChain()->getName();
      int resId2=a2->getResidue()->getId();
      string name2=a2->getName();
      string chainName2=a2->getResidue()->getChain()->getName();
      //use names to identify atoms in protein
      Atom *a3=molecule->getAtom(chainName1, resId1, name1);
      Atom *a4=molecule->getAtom(chainName2, resId2, name2);
      if( a3 == nullptr ) {
        cerr << "Cannot find atom with residue id " << resId1 << " and name " << name1 << endl;
        exit(-1);
      }
      if( a4 == nullptr ) {
        cerr << "Cannot find atom with residue id " << resId2 << " and name " << name2 << endl;
        exit(-1);
      }
//      makeCovBond(a3->getResidue(), a4->getResidue(), a3->getName(), a4->getName());
//      molecule->addCovBond(a3->getResidue(), a4->getResidue(), a3->getName(), a4->getName());
      molecule->addCovBond(a3, a4);
    }
  }

  //Create covalent bonds from CONECT records
  if( reference == nullptr ) {
    for( unsigned int i=0; i<conectRecords.size(); i++ ) {
      vector<int> row = conectRecords[i];
      Atom *a1=molecule->getAtom(row[0]);
      if( a1 == nullptr ) {
        cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with id " << row[0] << endl;
        exit(-1);
      }
      for(unsigned int j = 1; j < row.size(); j++ ){
        Atom *a2=molecule->getAtom(row[j]);
        if( a2 == nullptr ) {
          cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with id " << row[j] << endl;
          exit(-1);
        }
//        makeCovBond(a1->getResidue(), a2->getResidue(), a1->getName(), a2->getName());
//        molecule->addCovBond(a1->getResidue(), a2->getResidue(), a1->getName(), a2->getName());
        molecule->addCovBond(a1, a2);
//        cout << "Creating conect record bond between " << a1 << " and " << a2 << " in protein " << molecule->getName() << endl;
      }
    }
  }
  else {
    for( unsigned int i=0; i<conectRecords.size(); i++ ) {
      vector<int> row = conectRecords[i];
      Atom *a1 = reference->getAtom(row[0]);
      if (a1 == nullptr) {
        cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with id " << row[0] << endl;
        exit(-1);
      }
      int resId1 = a1->getResidue()->getId();
      string name1 = a1->getName();
      string chainName1 = a1->getResidue()->getChain()->getName();
      Atom *a3 = molecule->getAtom(chainName1, resId1, name1);
      if (a3 == nullptr) {
        cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with resi "<<resId1<<" and name "<<name1<<endl;
        exit(-1);
      }
      for (unsigned int j = 1; j < row.size(); j++) {
        Atom *a2 = reference->getAtom(row[j]);
        if (a2 == nullptr) {
          cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with id " << row[j] << endl;
          exit(-1);
        }
        int resId2 = a2->getResidue()->getId();
        string name2 = a2->getName();
        string chainName2 = a2->getResidue()->getChain()->getName();
        //use names to identify atoms in protein
        Atom *a4 = molecule->getAtom(chainName2, resId2, name2);
        if (a4 == nullptr) {
          cerr << "IO::readPdb - While parsing CONECT: Cannot find atom with resi "<<resId2<<" and name "<<name2<<endl;
          exit(-1);
        }
//        makeCovBond(a3->getResidue(), a4->getResidue(), a3->getName(), a4->getName());
//        molecule->addCovBond(a3->getResidue(), a4->getResidue(), a3->getName(), a4->getName());
        molecule->addCovBond(a3, a4);
        cout << "Creating conect record bond between " << a3 << " and " << a4 << " in protein " << molecule->getName()
             << endl;
      }
    }
  }

  // Fill in the second_cov_neighbor_list
  // Cannot do this step when the bond is still in creation because it won't know the neighbors of its neighbors yet
  for( auto const &atom: molecule->getAtoms()) {
    for( auto const &n1: atom->Cov_neighbor_list ) {
      for( auto const &n2: n1->Cov_neighbor_list ) {
        // check if n2 is ait itself. if yes, ignore it.
        if( n2 == atom )
          continue;
        // check whether n2 is already in the second_cov_neighbor_list of ait
        bool got_already=false;
        for( auto const &n3: atom->Second_cov_neighbor_list ) {
          if( n3 == n2 ) {
            got_already=true;
            break;
          }
        }
        if( !got_already )
          atom->Second_cov_neighbor_list.push_back(n2);
      }
    }

    //Print warning if atom has no covalent neighbors
    if( atom->Cov_neighbor_list.size() == 0 ) {
      cerr << "IO::readPdb - warning: Atom " << atom <<
      " has no covalent neighbors. Probably the atom-name isn't in the residue profiles" << endl;
      //exit(-1);
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
    if(atomName1<atomName2) tmp+=atomName1+"_"+atomName2;
    else                    tmp+=atomName2+"_"+atomName1;
    //log("debugRas")<<"Inserting "<<tmp<<endl;;
    fixedBonds.insert(tmp);
  }

  // Change the number of bars of locked bonds and peptide bonds to 6.
  // Cannot do this step when the bond is just created because it won't know the number of neighbors yet
  // No need to work on H-bonds because we assume all of them have 5 bars.
  for (auto const& bond: molecule->getCovBonds()) {

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

    // Set the number of bars to six for all these bonds!
    if(fixedBonds.find(query1)!=fixedBonds.end() || fixedBonds.find(query2)!=fixedBonds.end()){
      bond->Bars = 6;
    }

  }

//  // Assign main-chain/side-chain for all atoms
//  for (auto const& atom: molecule->getAtoms()){
//    string atom_name = atom->getName();
//    if (atom_name.compare("N")==0 || atom_name.compare("CA")==0 || atom_name.compare("C")==0 || atom_name.compare("O")==0) // this is a main-chain atom
//      atom->setAsMainchainAtom();
//    else if (
//        atom_name.compare("P")==0 || atom_name.compare("OP1")==0 || atom_name.compare("OP2")==0 ||
//        atom_name.compare("O5'")==0 ||
//        atom_name.compare("C5'")==0 ||
//        atom_name.compare("C4'")==0 || atom_name.compare("O4'")==0 ||
//        atom_name.compare("C3'")==0 ||
//        atom_name.compare("O3'")==0 ) // this is a main-chain atom in RNA
//      atom->setAsMainchainAtom();
//    else if (atom_name.find("H")!=string::npos) { // this is a H-atom
//      if(atom->Cov_neighbor_list.size() > 0) {
//        Atom* cov = atom->getIthCovNeighbor(0);
//        string cov_name = cov->getName();
//        if (cov_name.compare("N")==0 || cov_name.compare("CA")==0 || cov_name.compare("C")==0 || cov_name.compare("O")==0 ||
//            cov_name.compare("P")==0 || cov_name.compare("OP1")==0 || cov_name.compare("OP2")==0 ||
//            cov_name.compare("O5'")==0 ||
//            cov_name.compare("C5'")==0 ||
//            cov_name.compare("C4'")==0 ||	cov_name.compare("O4'")==0 ||
//            cov_name.compare("C3'")==0 ||
//            cov_name.compare("O3'")==0 )
//          atom->setAsMainchainAtom();
//      }
//      else
//        atom->setAsSidechainAtom();
//    }
//    else
//      atom->setAsSidechainAtom();
//  }


  //Create H-bonds from REMARK records
  for(const pair<int,int>& constraint: torsionConstraints) {
    Atom* hatom = molecule->getAtom(constraint.second);
    if (hatom == nullptr) {
      cerr << "IO::readPdb - Invalid atom-id specified: " << constraint.second << endl;
      exit(-1);
    }
    Atom* oatom = molecule->getAtom(constraint.first);
    if (oatom == nullptr) {
      cerr << "IO::readHbonds - Invalid atom-id specified: " << constraint.first << endl;
      exit(-1);
    }

    //Check if hatom and oatom were assigned correctly
    if (hatom->isHeavyAtom()) {
      oatom = hatom;
      hatom = molecule->getAtom(constraint.first);
    }

    //Assign donors and base atoms
    Atom* donor = hatom->getFirstCovNeighbor();
    Atom* AA = oatom->getFirstCovNeighbor();
    Hbond *new_hb = new Hbond(hatom, oatom, donor, AA, DEFAULT_HBOND_ENERGY);
    molecule->addHbond(new_hb);
  }

  //For backwards compatibility: If hbondFile set, read additional hydrogen bonds
  if(!hbondFile.empty()){
    readHbonds(hbondMethod, hbondFile, molecule);
  }

  molecule->sortHbonds();
  molecule->buildRigidBodies(movingResidues); //Necessary to do before building spanning tree
  molecule->buildSpanningTree(roots); //Necessary before conformations are defined
  molecule->setConfiguration(new Configuration(molecule));
  molecule->setCollisionFactor(1.0); //Sets the initial collisions
  return molecule;
}

//void IO::makeCovBond (Residue* res1, Residue* res2, string atom_name1, string atom_name2) {
//  //TODO: None of this needs to be in IO. Create a function in Molecule that does this and call it directly
//  Atom* atom1 = res1->getAtom(atom_name1);
//  Atom* atom2 = res2->getAtom(atom_name2);
//  bool valid = true;
//  if (atom1==nullptr) {
//    //cerr << "Warning: No such atom (" << atom_name1 << ") in Res" << res1->getId() << endl;
//    valid = false;
//  }
//  if (atom2==nullptr) {
//    //cerr << "Warning: No such atom (" << atom_name2 << ") in Res" << res2->getId() << endl;
//    valid = false;
//  }
//  if (valid) {
//    if(atom1->getId()>atom2->getId())
//      std::swap(atom1,atom2);
//
//    if(std::find(atom1->Cov_neighbor_list.begin(), atom1->Cov_neighbor_list.end(), atom2)==atom1->Cov_neighbor_list.end()) {
//      Bond *new_cb = new Bond(atom1, atom2, "COV");
//      res1->getChain()->getMolecule()->addCovBond(new_cb);
//    }
//  }
//}

void IO::readHbonds(const std::string& hbondMethod, const std::string& hbondFile, Molecule* mol)
{
  if(hbondMethod=="user")
    IO::readHbonds( mol, hbondFile );
  else if(hbondMethod=="rnaview")
    IO::readHbonds_rnaview( mol, hbondFile, true );
  else if(hbondMethod=="first" || hbondMethod=="FIRST")
    IO::readHbonds_first( mol, hbondFile );
  else if(hbondMethod=="kinari" || hbondMethod=="KINARI")
    IO::readHbonds_kinari( mol, hbondFile );
  else if(hbondMethod=="hbplus" || hbondMethod=="hbPlus")
    IO::readHbonds_hbPlus( mol, hbondFile );
  else if(hbondMethod=="vadar")
    IO::readHbonds_vadar( mol, hbondFile );
  else if(hbondMethod=="dssr")
    IO::readHbonds_dssr( mol, hbondFile );
  else if(hbondMethod=="identify")
    HbondIdentifier::identifyHbonds(mol);
}

//void IO::readDssp (Molecule * protein, string dssp_file) {
//  ifstream input(dssp_file.c_str());
//  if (!input.good()) {
//    cerr << "Error: cannot open file " << dssp_file << endl;
//    exit(1);
//  }
//  string sse_symbol="", last_symbol="x", chain_name = "", last_chain_name="";
//  int res_id, last_res_id=-100;
//  Chain* chain;
//  Residue* res;
//  string line;
//  bool start_read = false;
//  int helix_index=0, loop_index=0;
//  while ( getline(input,line) ) {
//    if (input.eof()) break;
//    if (line.at(2)=='#') {start_read = true; continue;}
//    if (!start_read) continue;
//    if (line.at(13)=='!') continue; // there is a discontinuity in either the backbone or the chain. skip this line
//
//    // Read in the symbol in DSSP file
//    chain_name = line.substr(11,1); // line[12]
//    if (chain_name.compare(last_chain_name)!=0) {
//      chain = protein->getChain(chain_name);
//      if (chain==nullptr) {
//        cerr << "Error: No such chain (" << chain_name << ") in m_molecule." << endl;
//        exit(1);
//      }
//      last_chain_name = chain_name;
//    }
//    res_id = atoi(line.substr(6,5).c_str()); // line[7:11]
//    res = chain->getResidue(res_id);
//    if (res==nullptr) {
//      cerr << "Error: No such residue " << res_id << " in chain (" << chain_name << ")." << endl;
//      exit(1);
//    }
//    sse_symbol = line.substr(16,1); // line[17]
//    if (sse_symbol.compare(" ")==0) // represents for loop, change the symbol to L
//      sse_symbol = "L";
//
//    // Convert the symbols and number the helices and loops
//    // Currently, no numbering on strands because no information on how to pair strands into sheets.
//    // G, H, I are all helices -> A.
//    // E and B are strands -> B.
//    // The rest are loops -> L.
//    if (sse_symbol.compare("G")==0 || sse_symbol.compare("H")==0 || sse_symbol.compare("I")==0) { // this is a helix
//      if (last_symbol.at(0)!='A' || last_res_id!=res_id-1)
//        ++helix_index;
//      sse_symbol = "A"+Util::i2s(helix_index);
//    }
//    else if (sse_symbol.compare("E")==0 || sse_symbol.compare("B")==0) { // this is a strand
//      sse_symbol = "B";
//    }
//    else {
//      if (last_symbol.at(0)!='L' || last_res_id!=res_id-1)
//        ++loop_index;
//      sse_symbol = "L"+Util::i2s(loop_index);
//    }
//
//    // Assign to the residue
//    res->setSSE(sse_symbol);
//    last_res_id = res_id;
//    last_symbol = sse_symbol;
//  }
//  input.close();
//}


//void IO::readRigidbody (Molecule * molecule) {
//  assert(molecule->getAtoms().size()>0);
//  //Create disjoint set
//  DisjointSets ds(molecule->getAtoms()[molecule->size() - 1]->getId() + 1); //Assumes the last atom has the highest id.
//
//  //For each atom, a1, with exactly one cov neighbor and not participating in an hbond, a2, call Union(a1,a2)
//  for (int i=0;i<molecule->size();i++){
//    Atom* atom = molecule->getAtoms()[i];
//    if(atom->Cov_neighbor_list.size()==1 && atom->Hbond_neighbor_list.size()==0){
//      ds.Union(atom->getId(), atom->Cov_neighbor_list[0]->getId());
//      //cout<<"Only one neighbor: "<<atom->getName()<<" "<<atom->getId()<<" - "<<atom->Cov_neighbor_list[0]->getName()<<" "<<atom->Cov_neighbor_list[0]->getId()<<endl;
//    }
//  }
//
////  //Create a map to determine if a covalent bond is fixed
////  set<string> fixedBonds;
////  for(int rp=0;;rp++){
////    if(strcmp(FIXED_BOND_PROFILES[rp][0],"END")==0) break;
////    string tmp(FIXED_BOND_PROFILES[rp][0]);
////    tmp+="_";
////    string atomName1(FIXED_BOND_PROFILES[rp][1]);
////    string atomName2(FIXED_BOND_PROFILES[rp][2]);
////    if(atomName1<atomName2)	tmp+=atomName1+"_"+atomName2;
////    else					tmp+=atomName2+"_"+atomName1;
////    //log("debugRas")<<"Inserting "<<tmp<<endl;;
////    fixedBonds.insert(tmp);
////  }
//
//  //For each fixed bond (a1,a2) call Union(a1,a2)
//  for (auto const& bond: molecule->getCovBonds()){
//
//    ///****************************************
//    ///UPDATED: we use two ways for now to identify locked bonds!
//    ///The isLocked/isPeptide bond check in Bond and the profiles
//    ///As for some profiles, it can be either locked or rotatable!
//
//    if( bond->Bars == 6){//This is fixed in the Bond -> isLocked
//      //cout<<"IO::readRigidBody - Bars=6"<<endl;
//      ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
//      continue;
//    }
//
////    Atom* a1 = bond->Atom1;
////    Atom* a2 = bond->Atom2;
////    //If the residue containing a1 (R1) precedes the one containing a2 (R2) there are four patterns in
////    //FIXED_BOND_PROFILES that can match: R1_a1_+a2, R1_+a2_a1, R2_-a1_a2 and R2_a2_-a1. Since atom names
////    //were ordered in the previous loop, however, we only have to check R1_+a2_a1 and R2_-a1_a2.
////    string query1, query2;
////    if(  a1->getResidue()->getNextResidue() == a2->getResidue()  ){
////      query1 = a1->getResidue()->getName()+"_+"+a2->getName()+"_"+a1->getName();
////      query2 = a2->getResidue()->getName()+"_-"+a1->getName()+"_"+a2->getName();
////    }else if(  a2->getResidue()->getNextResidue() == a1->getResidue()  ){
////      query1 = a1->getResidue()->getName()+"_+"+a1->getName()+"_"+a2->getName();
////      query2 = a2->getResidue()->getName()+"_-"+a2->getName()+"_"+a1->getName();
////    }else{
////      //If a1 and a2 are in the same residue simply put the lexicographically smallest name first.
////      if(a1->getName()<a2->getName())
////        query1 = a1->getResidue()->getName()+"_"+a1->getName()+"_"+a2->getName();
////      else
////        query1 = a1->getResidue()->getName()+"_"+a2->getName()+"_"+a1->getName();
////    }
////    //log("debugRas")<<bond<<" .. query: "<<query1<<" ("<<query2<<")"<<endl;
////
////    if(fixedBonds.find(query1)!=fixedBonds.end() || fixedBonds.find(query2)!=fixedBonds.end()){
////      ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
////      ///Set the number of bars to six for all these bonds!
////      bond->Bars = 6;
////    }
//  }
//
//
//  int c=0;
//  map<int,int> idMap;//Maps atom id's to rigid body id's for use in the DS structure.
//
//  //Map the set-ID's to RB-ID's and add bonded atoms to RBs.
//  for (int i=0;i<molecule->size();i++){
//    Atom* atom = molecule->getAtoms()[i];
//
//    //Map the set-id to the RB-id
//    int set_id = ds.FindSet(atom->getId());
//    int body_id;
//    if(idMap.find(set_id)!=idMap.end()) body_id = idMap.find(set_id)->second;
//    else {
//      body_id = c++;
//      idMap.insert( make_pair(set_id, body_id) );
//    }
//    //If the set containing a1 is not a rigid body: create one
//    if ( molecule->m_rigidBodyMap.find(body_id)==molecule->m_rigidBodyMap.end() ) {
//      Rigidbody* new_rb = new Rigidbody(body_id);
//      molecule->m_rigidBodyMap.insert( make_pair(body_id,new_rb) );
//    }
//    Rigidbody* rb = molecule->m_rigidBodyMap.find(body_id)->second;
//    if (!rb->containsAtom(atom)) rb->addAtom(atom);
//
//  }
//
//  //Delete small RBs and sort atoms within each RB
//  map<unsigned int, Rigidbody*>::iterator it = molecule->m_rigidBodyMap.begin();
//  while ( it!=molecule->m_rigidBodyMap.end() ) {
//    Rigidbody *rb = it->second;
//    if ( rb->Atoms.size() <= 0 ) {
//      cerr << "Error: rigid body " << rb->id() << " has no atoms." << endl;
//      molecule->m_rigidBodyMap.erase(it++);
//      delete rb;
//    }
//    else { // sort atom ids
//#ifndef WIN32 // Liangjun Zhang's tmp code
//      vector<Atom*>::iterator sit = rb->Atoms.begin();
//      vector<Atom*>::iterator eit = rb->Atoms.end();
//
//      sort(sit,eit,Atom::compare);
//#endif
//      // Determine if Rigidbody is on Mainchain
//      ++it;
//    }
//  }
//
//  //Store bonds in rigid bodies
//  for (auto const& bond: molecule->getCovBonds()){
//    int setId1 = ds.FindSet(bond->Atom1->getId());
//    molecule->m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
//    int setId2 = ds.FindSet(bond->Atom2->getId());
//    if(setId1!=setId2)
//      molecule->m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
//  }
//  for (auto const& bond: molecule->getHBonds()){
//    int setId1 = ds.FindSet(bond->Atom1->getId());
//    molecule->m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
//    int setId2 = ds.FindSet(bond->Atom2->getId());
//    if(setId1!=setId2)
//      molecule->m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
//  }
//
//
//}

//void IO::readRigidbody (Molecule * molecule, Selection& movingResidues) {
//  //Create disjoint set
//  DisjointSets ds(molecule->getAtoms()[molecule->size() - 1]->getId() + 1); //Assumes the last atom has the highest id.
//
//  //For each atom not in a residue in movingResidues call Union (rigidifies everything not in movingResidues)
//  for(auto const& chain: molecule->chains) {
//    Atom* lastAtom = nullptr;
//    for (auto const& res: chain->getResidues()) {
////      if(std::find(movingResidues.begin(), movingResidues.end(), res->getId())!=movingResidues.end())
//      if(movingResidues.inSelection(res)) {
//        log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Not rigidifying residue " << res->getId() << endl;
//        continue; //Skip residue if its in movingResidues
//      }else {
//        log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Rigidifying residue " << res->getId() << endl;
//      }
//
//      for (Atom *const &res_atom: res->getAtoms()){
//        if(lastAtom==nullptr) {
//          lastAtom = res_atom;
//          continue;
//        }
//        ds.Union(lastAtom->getId(), res_atom->getId());
//        log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Joining " << lastAtom->getId() << " - " << res_atom->getId() << endl;
//      }
//    }
//  }
//
//  //For each atom, a1, with exactly one cov neighbor and not participating in an hbond, a2, call Union(a1,a2)
//  for (int i=0;i<molecule->size();i++){
//    Atom* atom = molecule->getAtoms()[i];
//    if(atom->Cov_neighbor_list.size()==1 && atom->Hbond_neighbor_list.size()==0){
//      ds.Union(atom->getId(), atom->Cov_neighbor_list[0]->getId());
//      log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Joining " << atom->getId() << " - " << atom->Cov_neighbor_list[0]->getId() << endl;
//      //cout<<"Only one neighbor: "<<atom->getName()<<" "<<atom->getId()<<" - "<<atom->Cov_neighbor_list[0]->getName()<<" "<<atom->Cov_neighbor_list[0]->getId()<<endl;
//    }
//  }
//
////  //Create a map to determine if a covalent bond is fixed
////  set<string> fixedBonds;
////  for(int rp=0;;rp++){
////    if(strcmp(FIXED_BOND_PROFILES[rp][0],"END")==0) break;
////    string tmp(FIXED_BOND_PROFILES[rp][0]);
////    tmp+="_";
////    string atomName1(FIXED_BOND_PROFILES[rp][1]);
////    string atomName2(FIXED_BOND_PROFILES[rp][2]);
////    if(atomName1<atomName2)	tmp+=atomName1+"_"+atomName2;
////    else					tmp+=atomName2+"_"+atomName1;
////    //log("debugRas")<<"Inserting "<<tmp<<endl;;
////    fixedBonds.insert(tmp);
////  }
//
//  //For each fixed bond (a1,a2) call Union(a1,a2)
//  for (auto const& bond: molecule->getCovBonds()){
//
//    ///****************************************
//    ///UPDATED: we use two ways for now to identify locked bonds!
//    ///The isLocked/isPeptide bond check in Bond and the profiles
//    ///As for some profiles, it can be either locked or rotatable!
//
//    if( bond->Bars == 6){//This is fixed in the Bond -> isLocked
//      //cout<<"IO::readRigidBody - Bars=6"<<endl;
//      ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
//      log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Joining " << bond->Atom1->getId() << " - " << bond->Atom2->getId() << endl;
//      continue;
//    }
//
////    Atom* a1 = bond->Atom1;
////    Atom* a2 = bond->Atom2;
////    //If the residue containing a1 (R1) precedes the one containing a2 (R2) there are four patterns in
////    //FIXED_BOND_PROFILES that can match: R1_a1_+a2, R1_+a2_a1, R2_-a1_a2 and R2_a2_-a1. Since atom names
////    //were ordered in the previous loop, however, we only have to check R1_+a2_a1 and R2_-a1_a2.
////    string query1, query2;
////    if(  a1->getResidue()->getNextResidue() == a2->getResidue()  ){
////      query1 = a1->getResidue()->getName()+"_+"+a2->getName()+"_"+a1->getName();
////      query2 = a2->getResidue()->getName()+"_-"+a1->getName()+"_"+a2->getName();
////    }else if(  a2->getResidue()->getNextResidue() == a1->getResidue()  ){
////      query1 = a1->getResidue()->getName()+"_+"+a1->getName()+"_"+a2->getName();
////      query2 = a2->getResidue()->getName()+"_-"+a2->getName()+"_"+a1->getName();
////    }else{
////      //If a1 and a2 are in the same residue simply put the lexicographically smallest name first.
////      if(a1->getName()<a2->getName())
////        query1 = a1->getResidue()->getName()+"_"+a1->getName()+"_"+a2->getName();
////      else
////        query1 = a1->getResidue()->getName()+"_"+a2->getName()+"_"+a1->getName();
////    }
////    //log("debugRas")<<bond<<" .. query: "<<query1<<" ("<<query2<<")"<<endl;
////
////    if(fixedBonds.find(query1)!=fixedBonds.end() || fixedBonds.find(query2)!=fixedBonds.end()){
////      ds.Union(bond->Atom1->getId(), bond->Atom2->getId());
////      log("debug") << "IO::readRigidbody["<< __LINE__<<"] - Joining " << bond->Atom1->getId() << " - " << bond->Atom2->getId() << endl;
////      ///Set the number of bars to six for all these bonds!
////      bond->Bars = 6;
////    }
//  }
//
//
//  int c=0;
//  map<int,int> idMap;//Maps atom id's to rigid body id's for use in the DS structure.
//
//  //Map the set-ID's to RB-ID's and add bonded atoms to RBs.
//  for (int i=0;i<molecule->size();i++){
//    Atom* atom = molecule->getAtoms()[i];
//
//    //Map the set-id to the RB-id
//    int set_id = ds.FindSet(atom->getId());
//    int body_id;
//    if(idMap.find(set_id)!=idMap.end()) body_id = idMap.find(set_id)->second;
//    else {
//      body_id = c++;
//      idMap.insert( make_pair(set_id, body_id) );
//    }
//    //If the set containing a1 is not a rigid body: create one
//    if ( molecule->m_rigidBodyMap.find(body_id)==molecule->m_rigidBodyMap.end() ) {
//      Rigidbody* new_rb = new Rigidbody(body_id);
//      molecule->m_rigidBodyMap.insert( make_pair(body_id,new_rb) );
//    }
//    Rigidbody* rb = molecule->m_rigidBodyMap.find(body_id)->second;
//    if (!rb->containsAtom(atom)) rb->addAtom(atom);
//
//  }
//
//  //Delete small RBs and sort atoms within each RB
//  map<unsigned int, Rigidbody*>::iterator it = molecule->m_rigidBodyMap.begin();
//  while ( it!=molecule->m_rigidBodyMap.end() ) {
//    Rigidbody *rb = it->second;
//    if ( rb->Atoms.size() <= 0 ) {
//      cerr << "Error: rigid body " << rb->id() << " has no atoms." << endl;
//      molecule->m_rigidBodyMap.erase(it++);
//      delete rb;
//    }
//    else { // sort atom ids
//#ifndef WIN32 // Liangjun Zhang's tmp code
//      vector<Atom*>::iterator sit = rb->Atoms.begin();
//      vector<Atom*>::iterator eit = rb->Atoms.end();
//
//      sort(sit,eit,Atom::compare);
//#endif
//      // Determine if Rigidbody is on Mainchain
////      rb->setMainchainRb();
//      ++it;
//    }
//  }
//
//  //Store bonds in rigid bodies
////  for (list<Bond *>::iterator bit=molecule->m_covBonds.begin(); bit != molecule->m_covBonds.end(); ++bit) {
//  for (auto const& bond: molecule->getCovBonds()){
////    Bond* bond = *bit;
//    //    cout<<"IO::readRigidbody() - "<<bond->Atom1->getId()<<" "<<bond->Atom2->getId()<<endl;
//    int setId1 = ds.FindSet(bond->Atom1->getId());
//    molecule->m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
//    int setId2 = ds.FindSet(bond->Atom2->getId());
//    if(setId1!=setId2)
//      molecule->m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
//  }
////  for (list<Hbond *>::iterator bit=molecule->m_hBonds.begin(); bit != molecule->m_hBonds.end(); ++bit) {
//  for (auto const& bond: molecule->getHBonds()) {
//    int setId1 = ds.FindSet(bond->Atom1->getId());
//    molecule->m_rigidBodyMap.find( idMap.find(setId1)->second )->second->addBond(bond);
//    int setId2 = ds.FindSet(bond->Atom2->getId());
//    if(setId1!=setId2)
//      molecule->m_rigidBodyMap.find( idMap.find(setId2)->second )->second->addBond(bond);
//  }
//
//
//}

void IO::writePdb (Molecule * molecule, string output_file_name) {
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
  if(molecule->m_conf!=nullptr){
    Configuration* c = molecule->m_conf;
    output << "REMARK\tID = " << c->m_id << endl;
    if(c->getParent() != nullptr)
      output << "REMARK\tParent ID = " << c->getParent()->m_id << endl;
    output << "REMARK\tTree depth = " << c->m_treeDepth << endl;
    output<<"REMARK\tTree-path = ";
    int count=0;
    while(c->getParent() != nullptr && c->getParent() != c) { output << c->m_id << " "; c=c->getParent(); count++; }
    output << c->m_id << endl;
    //Use molecule->m_conf for all information
    output<<"REMARK\tDistance_initial = "<<setprecision(3)<<molecule->m_conf->m_distanceToIni<<endl;
    output<<"REMARK\tDistance to parent = "<<setprecision(3)<<molecule->m_conf->m_distanceToParent<<endl;
    output<<"REMARK\tDistance to target = "<<setprecision(6)<<molecule->m_conf->m_distanceToTarget<<endl;
    output<<"REMARK\tMax violation = "<<setprecision(3)<<molecule->m_conf->m_maxConstraintViolation<<endl;
    output<<"REMARK\tClash prevention = "<<setprecision(3)<<molecule->m_conf->m_usedClashPrevention<<endl;
    output<<"REMARK\tClash free dofs = "<<setprecision(3)<<molecule->m_conf->m_clashFreeDofs<<endl;
    output<<"REMARK\tOverall dofs = "<<setprecision(3)<<molecule->m_conf->getNumDOFs()<<endl;
    output<<"REMARK\tMin collision factor = "<<setprecision(3)<<molecule->m_conf->m_minCollisionFactor<<endl;
    output<<"REMARK\tVdw energy = "<<setprecision(6)<<molecule->m_conf->m_vdwEnergy<<endl;
  }
//  for (vector<Atom*>::iterator atom_itr=molecule->getAtoms().begin(); atom_itr != molecule->getAtoms().end(); ++atom_itr) {
  for(auto const& atom: molecule->getAtoms()){
//    Atom* atom = *atom_itr;
    Residue* res = atom->getResidue();
    char buffer[100];
//    int bigRBId = atom->getBiggerRigidbody()==nullptr?0:atom->getBiggerRigidbody()->id();
    sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s  ",
            atom->getId(),atom->getName().c_str(),
            res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
            atom->m_position.x,atom->m_position.y,atom->m_position.z,
            atom->getBFactor(),
//            bigRBId,
            atom->getType().c_str() );
    string line(buffer);
    output << line << endl;
  }
  output.close();
}

void IO::writeQ (Molecule *protein, Configuration* referenceConf, string output_file_name) {

  string myfile_s = output_file_name.substr(0,output_file_name.length()-4) + ".txt";
  ofstream myfile(myfile_s.c_str());
  if(!myfile.is_open()) {
    cerr<<"Cannot write to "<<myfile_s<<endl;
    exit(-1);
  }
  if(protein->m_conf->getNumDOFs() != referenceConf->getNumDOFs()){
    log("dominik")<<"Configurations don't have same dof number, not writing a q file."<<endl;
    return;
  }

  //Keep track of changes (in magnitude)
  for (auto const& edge: protein->m_spanningTree->Edges){
    if(edge->getBond()==nullptr) continue;

    int dof_id = edge->getDOF()->getIndex();
    int resId = edge->getBond()->Atom1->getResidue()->getId();
    int cycleDOF_id = edge->getDOF()->getCycleIndex();
    bool onBackbone=false;
    if( edge->getBond()->Atom1->isBackboneAtom() && edge->getBond()->Atom2->isBackboneAtom())
      onBackbone = true;
    int second = 1;
    if(cycleDOF_id == -1)
      second = 0;
    else if(protein->m_conf->getNullspace()->isCovBondRigid(cycleDOF_id) )
      second = 2;
    double absChangeI = formatRangeRadian(protein->m_conf->getGlobalTorsion(dof_id) -
                                              referenceConf->getGlobalTorsion(dof_id));
    myfile << dof_id <<" "<<resId <<" "<<second<<" "<<absChangeI <<" ";
    //myfile <<protein->m_conf->m_sumProjSteps[dof_id];
    myfile<<" "<<onBackbone<<endl;
  }
  myfile.close();
}

void IO::writeBondLengthsAndAngles (Molecule *molecule, string output_file_name) {
  string covBondFile = output_file_name + "_allCovBonds.txt";
  ofstream output(covBondFile.c_str());
//  for (list<Bond*>::iterator bond_itr=molecule->m_covBonds.begin(); bond_itr!=molecule->m_covBonds.end(); ++bond_itr) {
  for (auto const& bond: molecule->getCovBonds()){
//    Bond* bond = (*bond_itr);
    Math3D::Vector3 bondVec = bond->Atom1->m_position - bond->Atom2->m_position;
    output << bondVec.norm() << endl;
  }
  output.close();

  string edgeFile = output_file_name + "_allEdges.txt";
  ofstream output1(edgeFile.c_str());
  for (vector<KinEdge*>::iterator edge_itr=molecule->m_spanningTree->Edges.begin(); edge_itr!=molecule->m_spanningTree->Edges.end(); ++edge_itr) {
    Bond* bond = (*edge_itr)->getBond();
    if(bond==nullptr) continue;
    Math3D::Vector3 bondVec = bond->Atom1->m_position - bond->Atom2->m_position;
    output1 << bondVec.norm() << endl;
  }
  output1.close();

  string anchorEdgeFile = output_file_name + "_allAnchors.txt";
  ofstream output2(anchorEdgeFile.c_str());
  for (vector< pair<KinEdge*,KinVertex*> >::iterator edge_itr=molecule->m_spanningTree->m_cycleAnchorEdges.begin(); edge_itr!=molecule->m_spanningTree->m_cycleAnchorEdges.end(); ++edge_itr) {
    Bond* bond = edge_itr->first->getBond();
    Math3D::Vector3 bondVec = bond->Atom1->m_position - bond->Atom2->m_position;
    output2 << bondVec.norm() << endl;
  }
  output2.close();
}


void IO::writeCovBonds (Molecule *molecule, string output_file_name) {
  ofstream output(output_file_name.c_str());
//  for (list<Bond *>::iterator bond_itr=molecule->m_covBonds.begin(); bond_itr != molecule->m_covBonds.end(); ++bond_itr) {
  for (auto const& bond: molecule->getCovBonds()){
    output << right << setw(8) << bond->Atom1->getId();
    output << right << setw(8) << bond->Atom2->getId();
    output << right << setw(8) << bond->Bars << endl;
  }
  output.close();
}

//void IO::readCovBonds (Molecule *molecule, string in_file_name) {
//
//  ifstream cov(in_file_name.c_str());
//  if (!cov.good()) {
//    cerr << "Error: cannot open file " << in_file_name << endl;
//    exit(1);
//  }
//  string line, temp;
//  while (!cov.eof())
//  {
//    getline(cov,line);
//    // log("debugRas") << line;
//    if(line.size()== 0) // the last empty line
//      break;
//
//    string atom_name1 = line.substr(0,8);
//    string atom_name2 = line.substr(8,8);
//    Atom *atom1 = molecule->getAtom(atoi(atom_name1.c_str()));
//    Atom *atom2 = molecule->getAtom(atoi(atom_name2.c_str()));
//
//    Bond * new_cb = new Bond(atom1, atom2, "COV");
//    molecule->addCovBond(new_cb);
//  }
//  cov.close();
//}

void IO::writeHbondsIn (Molecule *molecule, string output_file_name) {
  ofstream output(output_file_name.c_str());

  int count=0;

  //Header line
//  output << "H-ID ";
//  output << "Acc-ID ";
//  output << "energy"<<endl;

  for (auto const& bond: molecule->getHBonds()){
    output << bond->Hatom->getId()<<" ";
    output << bond->Acceptor->getId()<<" ";
    output << bond->getIniEnergy()<<" "; //<<endl;
    output << bond->donorHybridization()<<" ";
    output << bond->acceptorHybridization()<<endl;
  }
  output.close();
}


void IO::writeHbonds (Molecule *molecule, string output_file_name) {
  ofstream output(output_file_name.c_str());

  int count=0;

  //Header line
  output << "H-bond_ID ";
  output << "H-chain ";
  output << "H-resi ";
  output << "H-resn ";
  output << "H-atomn ";
  output << "H-ID ";
  output << "Acc-chain ";
  output << "Acc-resi ";
  output << "Acc-resn ";
  output << "Acc-atomn ";
  output << "Acc-ID ";
  output << "energy ";
  output << "#bars ";
  output << "length ";
  output << "angle_D_H_A ";
  output << "angle_H_A_AA"<<endl;

  for (auto const& bond: molecule->getHBonds()){
    output << ++count<<" ";
    output << bond->Hatom->getResidue()->getChain()->getName()<<" ";
    output << bond->Hatom->getResidue()->getId()<<" ";
    output << bond->Hatom->getResidue()->getName()<<" ";
    output << bond->Hatom->getName()<<" ";
    output << bond->Hatom->getId()<<" ";
    output << bond->Acceptor->getResidue()->getChain()->getName()<<" ";
    output << bond->Acceptor->getResidue()->getId()<<" ";
    output << bond->Acceptor->getResidue()->getName()<<" ";
    output << bond->Acceptor->getName()<<" ";
    output << bond->Acceptor->getId()<<" ";
    output << bond->getIniEnergy()<<" ";
    output << bond->Bars<<" ";
    output << bond->getLength()<<" ";
    output << bond->getAngle_D_H_A()<<" ";
    output << bond->getAngle_H_A_AA()<<endl;
  }
  output.close();
}

void IO::writeHbondsChange (Configuration *conf, string output_file_name) {
  ofstream output(output_file_name.c_str());
  Molecule* molecule = conf->updatedMolecule();

  int count=0;
  //Header line
  output << "H-ID ";
  output << "Acc-ID ";
  output << "energy ";
  output << "#bars ";
  output << "length ";
  output << "angle_D_H_A ";
  output << "angle_H_A_AA ";
  output << "energy_change ";
  output << "length_change ";
  output << "D_H_A_change ";
  output << "H_A_AA_change ";
  output << "H-bond_ID" <<endl;

  for (auto const& bond: molecule->getHBonds()){
    double energy = bond->computeEnergy();
    output << bond->Hatom->getId()<<" ";
    output << bond->Acceptor->getId()<<" ";
    output << energy<<" ";
    int numBars = bond->Bars;
    output << numBars<<" ";
    output << bond->getLength()<<" ";
    output << bond->getAngle_D_H_A()<<" ";
    output << bond->getAngle_H_A_AA()<<" ";

    double energyChange = energy - bond->getIniEnergy();
    double distanceChange = bond->getLength() - bond->getIniLength();
    double H_A_AA_Change = formatRangeRadian( bond->getAngle_H_A_AA() - bond->getIniAngle_H_A_AA() );
    double D_H_A_Change = formatRangeRadian( bond->getAngle_D_H_A() - bond->getIniAngle_D_H_A() );

    output << energyChange<<" ";
    output << distanceChange<<" ";
    output << D_H_A_Change<<" ";
    output << H_A_AA_Change<<" ";
    output << ++count<<endl;
  }
  output.close();
}

void IO::readHbonds (Molecule *molecule, string hbond_file_name) {
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
//    if (input >> energy_s){
//      energy = atof(energy_s.c_str());
//    }else{
      energy = DEFAULT_HBOND_ENERGY;
//    }
    hatom_id = atoi(atom1_sid.c_str());
    oatom_id = atoi(atom2_sid.c_str());

    hatom = molecule->getAtom(hatom_id);
    if(hatom==nullptr){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<hatom_id<<endl; exit(-1); }
    oatom = molecule->getAtom(oatom_id);
    if(oatom==nullptr){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<oatom_id<<endl; exit(-1); }

    //Check if hatom and oatom were assigned correctly
    if( hatom->isHeavyAtom() ) {
      oatom = hatom;
      hatom = molecule->getAtom(oatom_id);
    }

    //Assign donors and base atoms
    donor = hatom->getFirstCovNeighbor();
    AA = oatom->getFirstCovNeighbor();
    Hbond * new_hb = new Hbond(hatom, oatom, donor, AA, energy);
    molecule->addHbond(new_hb);
  }
  input.close();
}

void IO::readAnnotations (Molecule *molecule, string annotation_file_name){
  ifstream input(annotation_file_name.c_str());
  if(!input.is_open()) return;

  Chain* chain = (molecule->chains[0]);
  int residues = chain->getResidues()[chain->getResidues().size()-1]->getId()+1;//Index of last residue + 1
  molecule->residueAnnotations = new int[residues];
  for(int i=0;i<residues;i++)
    molecule->residueAnnotations[i] = 0;

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
    molecule->residueAnnotations[residue] = annotation;
    //log("debugRas")<<"resAnn["<<residue<<"] = "<<annotation<<endl;
  }

}

void IO::readHbonds_dssr(Molecule * molecule, string dssrFile){
  cerr<<"IO::readHbonds_dssr(..) not yet implemented. Sorry"<<endl;
  exit(-1);
}

void IO::readHbonds_rnaview(Molecule * molecule, string file, bool fillAnnotations){
  int residues=0;
  if(fillAnnotations){
    for(auto const& chain: molecule->chains){
      int lastRes = chain->getResidues()[chain->getResidues().size()-1]->getId()+1;//Index of last residue + 1
      if(lastRes>residues) residues = lastRes;
    }
    molecule->residueAnnotations = new int[residues];
    for(int i=0;i<residues;i++)
      molecule->residueAnnotations[i] = 0;
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
//    cerr<<"IO::readHbonds_rnaview - TODO: check that chains are read correctly ("<<chain1<<", "<<chain2<<")"<<endl;
    int res1 = atoi(Util::trim(line.substr(13,6)).c_str());
    int res2 = atoi(Util::trim(line.substr(23,6)).c_str());

    if(fillAnnotations){
      if(res1>0 && res1<residues)	molecule->residueAnnotations[res1] = 1;
      if(res2>0 && res2<residues)	molecule->residueAnnotations[res2] = 1;
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
    Atom* donor = molecule->getAtom(chain1, res1, donorName);
    if(donor==nullptr) { cerr<<"IO::readHbonds_rnaview(..) Error: Couldnt find "<<res1<<"_"<<donorName<<endl;exit(-1);}
    Atom* acceptor = molecule->getAtom(chain2, res2, acceptorName);
    if(acceptor==nullptr) { cerr<<"IO::readHbonds_rnaview(..) Error: Couldnt find "<<res2<<"_"<<acceptorName<<endl;exit(-1);}
    Atom* AA = acceptor->getFirstCovNeighbor();
    Atom* hatom = nullptr;
    for(vector<Atom*>::iterator ait = donor->Cov_neighbor_list.begin(); ait != donor->Cov_neighbor_list.end(); ait++){
      Atom* a = *ait;
      if(a->m_element==atomH){ hatom = a; break; }
    }
    if(hatom==nullptr){//TODO: Manually create H-atom
      //int newId = *(m_molecule->atoms.end()-1)
      //Atom* atom = new Atom("H",,pos,occ,B);
      //m_molecule->addAtom(chain_name,res_name,res_id,atom);

      cerr<<"IO::readHbonds_rnaview(..). Error: Couldnt locate hydrogen in residue "<<res1<<endl;
      cerr<<"IO::readHbonds_rnaview(..). Error: "<<line<<endl;
      exit(-1);
    }
    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    molecule->addHbond(new_hb);
    log("so")<<"RNAView file: Created hydrogen bond between "<<hatom<<
      "("<<hatom->getId()<<") and "<<acceptor<<"("<<acceptor->getId()<<")"<<endl;
  }

  if(fillAnnotations){
    for(int i=0;molecule->residueAnnotations[i]==0;i++)
      molecule->residueAnnotations[i]=1;
    for(int i=residues-1;molecule->residueAnnotations[i]==0;i--)
      molecule->residueAnnotations[i]=1;
  }
  //enableLogger("debugRas");
  log("so")<<"ResidueAnnotations: ";
  for(int i=0;i<residues;i++){
    log("so")<<molecule->residueAnnotations[i];
  }
  log("so")<<endl;

}

void IO::readHbonds_first(Molecule * molecule, string file){
  ifstream input(file.c_str());
  if(!input.is_open()) {cerr<<"IO::readHbonds_first(..) Error: Couldnt open "<<file<<" for reading"<<endl; exit(-1);}

  string line;
  while( getline(input, line) ){
    if(line.at(line.length()-1)=='\r') line = line.substr(0,line.length()-1);

    int id1 = atoi(Util::trim(line.substr(0,8)).c_str());
    int id2 = atoi(Util::trim(line.substr(8,8)).c_str());

    Atom* hatom = molecule->getAtom(id1);
    Atom* acceptor = molecule->getAtom(id2);
    Atom* AA = acceptor->getFirstCovNeighbor();
    Atom* donor = hatom->getFirstCovNeighbor();

    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    molecule->addHbond(new_hb);
    log("verbose")<<"FIRST file: Created hydrogen bond between "<<hatom<<" and "<<acceptor<<endl;
  }
}

void IO::readHbonds_kinari (Molecule *molecule, string hbond_file_name) {
  ifstream input(hbond_file_name.c_str());
  std::string line;
  string atom1_sid, atom2_sid, energy_s;
  int hatom_id, oatom_id;
  double energy;
  string nullInput;
  Atom *hatom, *oatom, *donor, *AA;
  while(std::getline(input,line)){
    std::istringstream input(line);
    input >> nullInput;
    input >> atom1_sid;
    input >> nullInput;
    input >> atom2_sid;
    input >> energy_s;

//    if (input >> energy_s){
//      energy = atof(energy_s.c_str());
//    }else{
    energy = DEFAULT_HBOND_ENERGY;
//    }

    hatom_id = atoi(atom1_sid.c_str());
    oatom_id = atoi(atom2_sid.c_str());

    hatom = molecule->getAtom(hatom_id);
    if(hatom==nullptr){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<hatom_id<<endl; exit(-1); }
    oatom = molecule->getAtom(oatom_id);
    if(oatom==nullptr){ cerr<<"IO::readHbonds - Invalid atom-id specified: "<<oatom_id<<endl; exit(-1); }

    //Check if hatom and oatom were assigned correctly
    if( hatom->isHeavyAtom() ) {
      oatom = hatom;
      hatom = molecule->getAtom(oatom_id);
    }

    //Assign donors and base atoms
    donor = hatom->getFirstCovNeighbor();
    AA = oatom->getFirstCovNeighbor();
    Hbond * new_hb = new Hbond(hatom, oatom, donor, AA, energy);
    molecule->addHbond(new_hb);
  }
  input.close();
}

void IO::readHbonds_hbPlus(Molecule *protein, string hbond_file_name) {
  ifstream input(hbond_file_name.c_str());
  std::string line;
  string donorString,donorType, acceptorString, acceptorType, distanceHAString, hybridState, rest;
  int hatomID, acceptorID;
  Atom *hatom, *acceptor, *donor, *AA;
  int lineCount = 1;

  while(std::getline(input,line)) {
    if (lineCount < 9) {
      ++lineCount;
      continue;
    }
    std::istringstream input(line);
//    input >> donorString;
//    input >> donorType;
//    input >> acceptorString;
//    input >> acceptorType;
//    input >> rest;
//    input >> hybridState;
//    for (int i = 0; i < 3; i++) {
//      input >> rest;
//    }
//    input >> distanceHAString;

//    double distanceHA = atof(distanceHAString.c_str());
    double distanceHA = atof(line.substr(52,4).c_str());
    cout<<distanceHA<<endl;
    //Identify donor
    string donorChain = line.substr(0, 1);
    if (donorChain == "-") //Default from hb2, if chain is empty
//      donorChain = protein->chains[0]->getName();
      donorChain = " ";

    int donorResId = atoi( (line.substr(1,4)).c_str() );
    donorType = line.substr(9,4);
    donor = protein->getAtom(donorChain, donorResId,donorType );
    if(donor==NULL){ cerr<<"IO::readHbonds - Invalid donor specified: "<<donorString<<" "<<donorType<<endl; exit(-1); }

    //Identify acceptor
    string acceptorChain = line.substr(14,1);
    acceptorType = line.substr(23,4);
    if (acceptorChain == "-") //Default from hb2, if chain is empty
//      acceptorChain = protein->chains[0]->getName();
      acceptorChain = " ";

    int acceptorResId = atoi( (line.substr(15,4)).c_str() );
    acceptor = protein->getAtom(acceptorChain, acceptorResId,acceptorType );
    if(acceptor==NULL){ cerr<<"IO::readHbonds - Invalid acceptor specified: "<<acceptorString<<" "<<acceptorType<<endl; exit(-1); }

    //Identify hydrogen
    double minDistance = 99999;
    hatom = nullptr;
    for( auto neighbor : donor->Cov_neighbor_list ){
      if( !neighbor->isHeavyAtom()){
        double distance = acceptor->distanceTo(neighbor);
        distance = distance*distance - distanceHA*distanceHA;
        if( distance < minDistance ){
          minDistance = distance;
          hatom = neighbor;
        }
      }
    }

    if(hatom==nullptr){ cerr<<"IO::readHbonds - Could not find suitable hydrogen for given donor and acceptor: "<<donorString<<" "<<donorType<<" ; "<<acceptorString<<" "<<acceptorType<<" Skipping."<<endl; continue; }

    AA = acceptor->getFirstCovNeighbor();
    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    protein->addHbond(new_hb);

    lineCount++;
  }
  input.close();
}


void IO::readHbonds_vadar(Molecule * molecule, string file){
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

    Atom* donor = molecule->getAtom(chain1, resId1, name1);
    if(donor==nullptr) { cerr<<"Cannot find res "<<resId1<<"/"<<name1<<endl; exit(-1); }
    Atom* acceptor = molecule->getAtom(chain1, resId2, name2);
    if(acceptor==nullptr) { cerr<<"Cannot find res "<<resId2<<"/"<<name2<<endl; exit(-1); }
    cout<<"Hbond "<<donor->getId()<<" "<<acceptor->getId()<<endl;
    Atom* hatom = donor->getFirstCovNeighbor();
    Atom* AA = acceptor->getFirstCovNeighbor();
    //Atom* hatom = donor->getFirstCovHydrogenNeighbor();
    //Atom* AA = acceptor->getFirstCovNeighbor();
    //if(hatom==nullptr){
    //	hatom = donor->getFirstCovHydrogenNeighbor();
    //	h
    //}
    //Atom* AA = hatom->getFirstCovNeighbor();

    Hbond * new_hb = new Hbond(hatom, acceptor, donor, AA);
    molecule->addHbond(new_hb);
    log("verbose")<<"vadar file: Created hydrogen bond between "<<hatom<<" and "<<acceptor<<endl;
  }

}

void IO::writePyMolScript(Molecule * rigidified, string pdb_file, string output_file_name, Molecule* iniMolecule) {

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
//  map< unsigned int, string > mapClusterIDtoColor;
  //map<unsigned int,Rigidbody*>::iterator rbit = molecule->m_conf->m_biggerRBMap.begin();
  //vector< pair< int, unsigned int> >::iterator sit = molecule->m_conf->m_sortedRBs.begin();
  Color::next_rcd_color_as_name(true);

//  while(sit != molecule->m_conf->m_sortedRBs.end() ){
//    if( sit->first >= MIN_CLUSTER_SIZE ){
  for(auto const& rb : rigidified->getRigidbodies()){
    if(rb->Atoms.size()>MIN_CLUSTER_SIZE){

      string color = Color::next_rcd_color_as_name();
      //mapClusterIDtoColor[sit->second] = color;
      //pymol_script << "color " << color << ", ( b > " << float(sit->second-0.01)
      //             << " and b < " << float(sit->second+0.01) << ")" << endl;
//      mapClusterIDtoColor[rb->id()] = color;
      pymol_script << "color " << color << ", ( b > " << float(rb->id()-0.01)
                   << " and b < " << float(rb->id()+0.01) << ")" << endl;
    }
//    sit++;
  }

  // Python commands to draw hbonds
  //////////////////////////////////////////////////////////////////////
  pymol_script << "# Draw hbonds as distance objects" << endl;
  pymol_script << "set dash_gap, 0.1" << endl;

  int site_1, site_2;
  vector< pair<KinEdge*,KinVertex*> >::iterator eit;

  if(iniMolecule) {
    for (eit = iniMolecule->m_spanningTree->m_cycleAnchorEdges.begin();
         eit != iniMolecule->m_spanningTree->m_cycleAnchorEdges.end(); eit++) {

      site_1 = eit->first->getBond()->Atom1->getId();
      site_2 = eit->first->getBond()->Atom2->getId();
      pymol_script << "distance allHbonds = id " << site_1 << " , id " << site_2 << endl;

    }
    pymol_script << "color red, allHbonds" << endl;
    pymol_script << "hide labels, allHbonds" << endl;
  }

  for (eit=rigidified->m_spanningTree->m_cycleAnchorEdges.begin(); eit != rigidified->m_spanningTree->m_cycleAnchorEdges.end(); eit++) {

    site_1 = eit->first->getBond()->Atom1->getId();
    site_2 = eit->first->getBond()->Atom2->getId();
    pymol_script << "distance hbondConstraints = id " << site_1 << " , id " << site_2 << endl;

  }
  pymol_script << "color yellow, hbondConstraints" << endl;
  pymol_script << "hide labels, hbondConstraints" << endl;

  // Create the rigid cluster objects for pymol, and color them. Only those
  // clusters larger than the min_output_cluster_size will have objects
  // created for them. The remaining clusters will be binned into bulk
  // objects.
  //////////////////////////////////////////////////////////////////////

  render_style = "lines";

  // for each rigid cluster, create a new pymol object of the cluster.
  unsigned int total_RC_objects = 0;
//  rbit = molecule->m_conf->m_biggerRBMap.begin();
//  sit = molecule->m_conf->m_sortedRBs.begin();
//  while( sit != molecule->m_conf->m_sortedRBs.end() ){
//    if( (sit->first) >= MIN_CLUSTER_SIZE ){
//      pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << sit->second<< " has "
//        << sit->first << " atoms." << endl;
//      pymol_script << "select RC" << sit->second << ", ( b > " << float(sit->second-0.01)
//        << " and b < " << float(sit->second+0.01) << ")" << endl;
//      pymol_script << "show " << render_style << ", RC" << sit->second << endl;
//      pymol_script << "color " << mapClusterIDtoColor[sit->second] << ", RC" << sit->second << endl << endl;
//    }
//    sit++;
//  }

//  for(auto const& it : molecule->m_spanningTree->Vertex_map) {
//    Rigidbody *rb = it.second->m_rigidbody;
//    if (rb != nullptr && rb->Atoms.size() > MIN_CLUSTER_SIZE) {
//      pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << rb->id()<< " has "
//                   << rb->Atoms.size() << " atoms." << endl;
//      pymol_script << "select RC" << rb->id() << ", ( b > " << float(rb->id()-0.01)
//                   << " and b < " << float(rb->id()+0.01) << ")" << endl;
//      pymol_script << "show " << render_style << ", RC" << rb->id() << endl;
//      pymol_script << "color " << mapClusterIDtoColor[rb->id()] << ", RC" << rb->id() << endl << endl;
//
//    }
//  }

  //Create biggest cluster as separate object
  pymol_script << "create biggestCluster, b < 0.99" << endl;


  // Some final global attributes to set.
  //////////////////////////////////////////////////////////////////////
  pymol_script << "bg_color white" << endl;
//  pymol_script << "clip slab, 200" << endl;
  pymol_script << "orient" << endl;

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
  if(protein->m_conf!=nullptr){
    Configuration* c = protein->m_conf;
    for(auto const& it: c->getMolecule()->m_spanningTree->Vertex_map){
      Rigidbody* rb = it.second->m_rigidbody;
      if(rb==nullptr) continue;
      for(auto const& atom: rb->Atoms){
        output<<atom->getId()<<endl;
      }
      output<<endl<<"NaN"<<endl<<endl;
    }
//    vector< pair<int, unsigned int> >::iterator it=c->m_sortedRBs.begin();
//
//    while( it != c->m_sortedRBs.end() ){
//
//      vector<Atom*>::iterator ait = c->m_biggerRBMap[it->second]->Atoms.begin();
//      while( ait != c->m_biggerRBMap[it->second]->Atoms.end() ){
//      output << (*ait)->getId() << endl;
//      ait++;
//    }
//      output << endl;
//      output << "NaN"<<endl;
//      output << endl;
//      it++;
//    }

  }
  output.close();

}


void IO::writeStats(Molecule * protein, string output_file_name, Molecule* rigidified){

  ofstream output( output_file_name.c_str() );
  if(!output.is_open()) {
    cerr<<"Cannot write to "<<output_file_name<<". You might need to create output directory first"<<endl;
    exit(-1);
  }
  if(protein->m_conf!=nullptr) {
    int diff = protein->m_spanningTree->getNumDOFs() - protein->m_spanningTree->getNumCycleDOFs();
    int sum = protein->m_conf->getNullspace()->getNullspaceSize() + diff;
    //int sum = m_molecule->m_conf->CycleNullSpace->getNullspace()Size + diff;

    output << "************* Statistics on the molecular structure *************" << endl;
    output << "Number of chains: " << protein->chains.size() << endl;
    output << "Number of atoms: " << protein->getAtoms().size() << endl;
    output << "Number of covalent bonds: " << protein->getCovBonds().size() << endl;
    output << "Number of hydrogen bonds: " << protein->m_spanningTree->m_cycleAnchorEdges.size() << endl;
    output << "Number of dihedrals in spanning tree: " << protein->m_spanningTree->getNumDOFs() << endl;
    output << "Number of free DOFs: " << diff << endl;
    output << "Number of cycle DOFs: " << protein->m_spanningTree->getNumCycleDOFs() << endl << endl;
    output << "************* Statistics on rigidity analysis *************" << endl;
    output << "Number of internal m_dofs (nullspace dimension): " << protein->m_conf->getNullspace()->getNullspaceSize()
           << endl;
    output << "Overall number of m_dofs (free + coordinated): " << sum << endl;
    output << "Number of rigidified covalent bonds: " << protein->m_conf->getNullspace()->getNumRigidDihedrals()
           << endl;
    output << "Number of rigidified hydrogen bonds: " << protein->m_conf->getNullspace()->getNumRigidHBonds() << endl;

    if (rigidified) { ///rigidified is a protein with collapsed rigid bodies
      output << "Number of rigid clusters: " << rigidified->m_conf->m_numClusters << endl;
      output << "Size of biggest cluster: " << rigidified->m_conf->m_maxSize << endl;
      output << "Index of biggest cluster: " << rigidified->m_conf->m_maxIndex << endl;
    } else { ///otherwise use normal protein
      output << "Number of rigid clusters: " << protein->m_conf->m_numClusters << endl;
      output << "Size of biggest cluster: " << protein->m_conf->m_maxSize << endl;
      output << "Index of biggest cluster: " << protein->m_conf->m_maxIndex << endl;
    }
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

  if(molecule->m_conf!=nullptr){

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
    //		molecule->setConfiguration(c);
    c->updateMolecule();
    //		molecule->alignToReference()
    c = molecule->m_conf;
    currId = c->m_id;
    Configuration* iniConf = c;
    output << "REMARK\tInitial Distance to target was "<<setprecision(6)<<c->m_distanceToTarget<<endl;

    molecule->writeRigidbodyIDToBFactor();

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
      //			for( eit = molecule->m_spanningTree->Edges.begin(); eit != molecule->m_spanningTree->Edges.end(); eit++){
      //				KinEdge* e = (*eit);
      //				int dofId = e->DOF_id;
      //				CTKResidue* res = e->Bond->Atom1->m_parentResidue;
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

      for (auto const& atom: molecule->getAtoms()) {
        Residue* res = atom->getResidue();
        char buffer[100];
        sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6.1f          %2s  ",
            atom->getId(),atom->getName().c_str(),
            res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
            //atom->Position.x,atom->Position.y,atom->Position.z,atom->getType().c_str());
          atom->m_position.x,atom->m_position.y,atom->m_position.z, atom->getBFactor(),atom->getType().c_str() );
//        atom->m_position.x,atom->m_position.y,atom->m_position.z,atom->getBiggerRigidbody()->id(),atom->getType().c_str() );
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

        for(auto const& cit: c->getChildren()){
          if(cit->m_id != path[i+1])
            continue;
          else{
            //						molecule->setConfiguration(*cit);
            cit->updateMolecule();
            break;
          }
        }
        c = molecule->m_conf;
        currId = c->m_id;
      }
    }
    delete[] path;
  }
  if( target != nullptr && target->m_conf != nullptr){//also add the additional path of the reversePlanner

    Configuration* c = target->m_conf;
    int currId = c->m_id;
    int treeDepth = c->m_treeDepth;

    Selection allSel;
    metrics::RMSD rmsd(allSel);
    double alignVal = rmsd.align(molecule,target);
    double distance = rmsd.distance_noOptimization(cFinalFwd,c);

    output << "REMARK\tAdding reversePlanning path of "<<treeDepth<<" samples."<<endl;
    output << "REMARK\tOverall path length is "<<firstPathLength+treeDepth<<" samples."<<endl;
    output << "REMARK\tStarting from reverse sample "<<currId<<endl;
    output << "REMARK\tRMSD at transition is "<<setprecision(3)<<distance<<endl;
    output << "REMARK\tReverse path follows "<<currId;

    while( c->getParent() != nullptr){
      c=c->getParent();
      currId = c->m_id;
      output << ", "<<currId;
    }
    output <<endl;
    c = target->m_conf;
    currId = c->m_id;
    Configuration* iniConf = c;

//    //Random permutation of all rigidbody ids for use in coloring
//    vector<int> rbidPerm;
//    for(size_t i=0;i<molecule->m_spanningTree->Vertex_map.size();i++)
//      rbidPerm.push_back(i);
//    std::random_shuffle ( rbidPerm.begin(), rbidPerm.end() );

    target->writeRigidbodyIDToBFactor();

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
      //			for( eit = target->m_spanningTree->Edges.begin(); eit != target->m_spanningTree->Edges.end(); eit++){
      //				KinEdge* e = (*eit);
      //				int dofId = e->DOF_id;
      //				CTKResidue* res = e->Bond->Atom1->m_parentResidue;
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

      for (auto const& atom: target->getAtoms()) {
        Residue* res = atom->getResidue();
        char buffer[100];
        //sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  ",
        sprintf(buffer,"ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s  ",
            atom->getId(),atom->getName().c_str(),
            res->getName().c_str(),res->getChain()->getName().c_str(),res->getId(),
            //atom->Position.x,atom->Position.y,atom->Position.z,atom->getType().c_str());
          atom->m_position.x,atom->m_position.y,atom->m_position.z,atom->getBFactor(),atom->getType().c_str() );
//        atom->m_position.x,atom->m_position.y,atom->m_position.z,atom->getBiggerRigidbody()->id(),atom->getType().c_str() );
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
      if( c->getParent() != nullptr){
        c=c->getParent();
        //				target->setConfiguration(c);
        c->updateMolecule();
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
//  map<unsigned int,Rigidbody*>::iterator rbit = molecule->m_conf->m_biggerRBMap.begin();
//  vector< pair< int, unsigned int> >::iterator sit = molecule->m_conf->m_sortedRBs.begin();
  Color::next_rcd_color_as_name(true);

//  while(sit != molecule->m_conf->m_sortedRBs.end() ){
//    if( sit->first >= MIN_CLUSTER_SIZE ){
//      string color = Color::next_rcd_color_as_name();
//      mapClusterIDtoColor[sit->second] = color;
//      pymol_script << "color " << color << ", ( b > " << float(sit->second-0.01)
//        << " and b < " << float(sit->second+0.01) << ")" << endl;
//    }
//    sit++;
//  }
  for(auto const& rb: molecule->getRigidbodies()){
    if(rb->Atoms.size()<MIN_CLUSTER_SIZE) continue;
    string color = Color::next_rcd_color_as_name();
    mapClusterIDtoColor[rb->id()] = color;
    pymol_script << "color " << color << ", ( b > " << float(rb->id()-0.01)
                 << " and b < " << float(rb->id()+0.01) << ")" << endl;
  }

  // Python commands to draw hbonds
  //////////////////////////////////////////////////////////////////////
  pymol_script << "# Draw hbonds as distance objects" << endl;
  pymol_script << "set dash_gap, 0.1" << endl;

  int site_1, site_2;
  vector< pair<KinEdge*,KinVertex*> >::iterator eit;

  for (auto  const& eit: molecule->m_spanningTree->m_cycleAnchorEdges){

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
//  rbit = molecule->m_conf->m_biggerRBMap.begin();
//  sit = molecule->m_conf->m_sortedRBs.begin();
//  while( sit != molecule->m_conf->m_sortedRBs.end() ){
//    if( (sit->first) >= MIN_CLUSTER_SIZE ){
//  pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << sit->second<< " has "
//               << sit->first << " atoms." << endl;
//  pymol_script << "select RC" << sit->second << ", ( b > " << float(sit->second-0.01)
//               << " and b < " << float(sit->second+0.01) << ")" << endl;
//  pymol_script << "show " << render_style << ", RC" << sit->second << endl;
//  pymol_script << "set line_width = 3, " << "RC" << sit->second << endl;
//  pymol_script << "color " << mapClusterIDtoColor[sit->second] << ", RC" << sit->second << endl << endl;
//}
//sit++;
//}
  for(auto const& rb: molecule->getRigidbodies()){
    if(rb->Atoms.size()>=MIN_CLUSTER_SIZE){
      pymol_script << "# Rigid Cluster # " << ++total_RC_objects << " and ID " << rb->id()<< " has "
        << rb->Atoms.size() << " atoms." << endl;
      pymol_script << "select RC" << rb->id() << ", ( b > " << float(rb->id()-0.01)
        << " and b < " << float(rb->id()+0.01) << ")" << endl;
      pymol_script << "show " << render_style << ", RC" << rb->id() << endl;
      pymol_script << "set line_width = 3, " << "RC" << rb->id() << endl;
      pymol_script << "color " << mapClusterIDtoColor[rb->id()] << ", RC" << rb->id() << endl << endl;
    }
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



std::vector< std::tuple<Atom*, Atom*, double> > IO::readRelativeDistances(const std::string& fname, Molecule* mol){
  ifstream input(fname.c_str());
  if(!input.is_open()) {cerr<<"IO::readrelativeDistances(..) - Error: Couldnt open "<<fname<<" for reading"<<endl; exit(-1);}

  std::vector< std::tuple<Atom*, Atom*, double> > ret;
  string line;
  while( getline(input,line) ){
    vector<string> tokens = Util::split(line, ",");
    if(tokens.size()!=3){
      cerr<<"IO::readRelativeDistances("<<fname<<") - error: Each restraint file line must have 2 atom selections and 1 distance that are comma-separated"<<endl;
      exit(-1);
    }

    Selection s1(Util::trim(tokens[0]));
    Selection s2(Util::trim(tokens[1]));
    double dist = std::stod(Util::trim(tokens[2]));

    Atom* a1 = nullptr;
    if(s1.getSelectedResidues(mol).size()==1) {
      a1 = s1.getSelectedResidues(mol)[0]->getAtom("CA");
      if(a1==nullptr) a1 = s1.getSelectedResidues(mol)[0]->getAtoms().front();
    }else{
      if(s1.getSelectedAtoms(mol).size()==1) {
        a1 = s1.getSelectedAtoms(mol).at(0);
      }else{
        cerr<<"IO::readRelativeDistances("<<fname<<") - error: Restraint selections must match exactly 1 atom or 1 residue"<<endl;
        exit(-1);
      }
    }

    Atom* a2 = nullptr;
    if(s2.getSelectedResidues(mol).size()==1) {
      a2 = s2.getSelectedResidues(mol)[0]->getAtom("CA");
      if(a2==nullptr) a2 = s2.getSelectedResidues(mol)[0]->getAtoms().front();
    }else{
      if(s2.getSelectedAtoms(mol).size()==1) {
        a2 = s2.getSelectedAtoms(mol).at(0);
      }else{
        cerr<<"IO::readRelativeDistances("<<fname<<") - error: Restraint selections must match exactly 1 atom or 1 residue"<<endl;
        exit(-1);
      }
    }


    ret.push_back( std::make_tuple( a1,a2, dist ) );
  }

  return ret;
}

void IO::writeNewSample(Configuration *conf, Configuration *ref, int sample_num, const string &workingDir, int saveData) {
//	const string& out_path = ExploreOptions::getOptions()->workingDirectory;
  const string &out_path = workingDir;
  const string &name = conf->getMolecule()->getName();
  string out_file = out_path + "output/" + name + "_new_" + std::to_string(static_cast<long long>(sample_num)) + ".pdb";

  if (saveData > 0) {

    Molecule *protein = conf->updatedMolecule();
    IO::writePdb(protein, out_file);
  }


  if (saveData > 1) {
    Molecule *protein = conf->updatedMolecule();

    string out_q = out_path + "output/" + name + "_q_" + std::to_string(static_cast<long long>(sample_num)) + ".txt";

    IO::writeQ(protein, ref, out_q);
  }


  if (saveData > 2) {
    Molecule *protein = conf->updatedMolecule();

    // Save Jacobian and Nullspace to file
    string outJac = out_path + "output/" + name + "_jac_" +
                    std::to_string(static_cast<long long>(sample_num))
                    + ".txt";
    string outNull = out_path + "output/" + name + "_nullSpace_" +
                     std::to_string(static_cast<long long>(sample_num))
                     + ".txt";
    // Save singular values
    string outSing = out_path + "output/" + name + "_singVals_" +
                     std::to_string(static_cast<long long>(sample_num))
                     + ".txt";
    // Save pyMol coloring script
    string pyMol = out_path + "output/" + name + "_pyMol_" +
                   std::to_string(static_cast<long long>(sample_num))
                   + ".pml";
    string rbFile = out_path + "output/" + name + "_RBs_" +
                    std::to_string(static_cast<long long>(sample_num))
                    + ".txt";
    string covFile = out_path + "output/" + name + "_covBonds_" +
                     std::to_string(static_cast<long long>(sample_num))
                     + ".txt";
    string statFile = out_path + "output/" + name + "_stats_" +
                      std::to_string(static_cast<long long>(sample_num))
                      + ".txt";

    IO::writePyMolScript(protein, out_file, pyMol);

    ///Write Jacobian
    //gsl_matrix_outtofile(conf->CycleJacobian, outJac);
    ///Write Null Space matrix
    //gsl_matrix_outtofile(conf->CycleNullSpace->V,outNull);
    //gsl_vector_outtofile(conf->CycleNullSpace->singularValues,outSing);

    if (NullspaceSVD *derived = dynamic_cast<NullspaceSVD *>(conf->getNullspace())) {
      derived->writeMatricesToFiles(outJac, outNull, outSing);
    }
    IO::writeRBs(protein, rbFile);
    IO::writeStats(protein, statFile);

  }

}
