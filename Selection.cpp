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
#include <regex>
#include "Selection.h"
#include "Logger.h"
#include "Util.h"

using namespace std;


Selection::Selection() : Selection("all") 
{
}

Selection::Selection(const string& selectionPattern) :
    m_selectionPattern(selectionPattern),
    m_rootClause(parseClause(selectionPattern))
{
}

vector<Atom*>& Selection::getSelectedAtoms( const Molecule* mol )
{
  if( m_cachedAtoms.count(mol)==0 ){
    m_cachedAtoms[mol] = std::move(std::vector<Atom*>());
    vector<Atom*>& ret = m_cachedAtoms[mol];

    for(auto const& a: mol->getAtoms()){
      if( m_rootClause->inSelection(a) )
        ret.push_back(a);
    }
  }

  return m_cachedAtoms[mol];
}

vector<Residue*>& Selection::getSelectedResidues( const Molecule *mol )
{
  if(m_cachedResidues.count(mol)==0) {

    m_cachedResidues[mol] = std::move(std::vector<Residue*>());
    vector<Residue *>& ret = m_cachedResidues[mol];
    for (auto chain: mol->chains) {
      for (auto res: chain->getResidues()) {

        //Test whether all atoms in res are selected
        bool includeRes = true;
        for (auto atom: res->getAtoms()) {
          if (!m_rootClause->inSelection(atom)) {
            includeRes = false;
            break;
          }
        }

        if (includeRes) ret.push_back(res);
      }
    }
  }

  return m_cachedResidues[mol];
}

std::vector<Bond *>& Selection::getSelectedBonds( const Molecule *mol )
{
  if(m_cachedBonds.count(mol)==0) {
    m_cachedBonds[mol] = std::move(vector<Bond*>());
    vector<Bond *>& ret = m_cachedBonds[mol];
    for (auto const &bond: mol->getCovBonds()) {
      //Check if both end-atoms are in selection
      if (m_rootClause->inSelection(bond->Atom1) && m_rootClause->inSelection(bond->Atom2))
        ret.push_back(bond);
    }
  }

  return m_cachedBonds[mol];
}

bool Selection::inSelection( const Atom* a ) const { return m_rootClause->inSelection(a); }

bool Selection::inSelection( const Bond* b ) const{
  return m_rootClause->inSelection(b->Atom1) && m_rootClause->inSelection(b->Atom2);
}

bool Selection::inSelection( const Residue* r) const{
  for(auto const& a: r->getAtoms())
    if(!m_rootClause->inSelection(a)) return false;
  return true;
}


Selection::Clause* Selection::parseClause(const std::string &input) {
  if(Util::contains(input, " or ") || Util::contains(input, " + ")) return new OrClause(input);
  if(Util::contains(input, " and ") || Util::contains(input, " & ")) return new AndClause(input);
  if(Util::startsWith(input, "not ") || Util::startsWith(input, "!")) return new NotClause(input);
  if(Util::startsWith(input, "resi ")) return new ResiClause(input);
  if(Util::startsWith(input, "resn ")) return new ResnClause(input);
  if(Util::startsWith(input, "name ")) return new NameClause(input);
  if(Util::startsWith(input, "elem ")) return new ElemClause(input);
  if(input=="all")      return new AllClause(input);
  if(input=="")         return new AllClause(input);
  if(input=="heavy")    return new HeavyClause(input);
  if(input=="hydro")    return new HydroClause(input);
  if(input=="backbone") return new BackboneClause(input);

  throw std::runtime_error("KGS error: Unrecognized pattern: "+input);
}

Selection::OrClause::OrClause(const std::string& input)
{
  vector<string> subClausesPlus = Util::split(input, " + ");
  for(const string& subClausePlus: subClausesPlus){
    vector<string> subClausesOr = Util::split(subClausePlus, " or ");
    for(const string& subClauseOr: subClausesOr){
      m_childClauses.push_back( parseClause(subClauseOr) );
    }
  }
}
bool Selection::OrClause::inSelection(const Atom* a) const
{
  for(const Clause* subClause: m_childClauses)
    if( subClause->inSelection(a) ) return true;
  return false;
}

Selection::AndClause::AndClause(const std::string& input)
{
  vector<string> subClausesAmp = Util::split(input, " & ");
  for(const string& subClauseAmp: subClausesAmp){
    vector<string> subClausesAnd = Util::split(subClauseAmp, " and ");
    for(const string& subClauseAnd: subClausesAnd){
      m_childClauses.push_back( parseClause(subClauseAnd) );
    }
  }
}
bool Selection::AndClause::inSelection(const Atom* a) const
{
  for(const Clause* subClause: m_childClauses)
    if( !subClause->inSelection(a) ) return false;
  return true;
}

Selection::NotClause::NotClause(const std::string& input)
{
  if(Util::startsWith(input, "not ")) m_childClause = parseClause(input.substr(4));
  else if(Util::startsWith(input, "! ")) m_childClause = parseClause(input.substr(2));
  else if(Util::startsWith(input, "!")) m_childClause = parseClause(input.substr(1));
}
bool Selection::NotClause::inSelection(const Atom* a) const
{
  return !(m_childClause->inSelection(a));
}

Selection::ResiClause::ResiClause(const std::string& input)
{
  regex intListPat("\\-?[[:digit:]]+(\\+\\-?[[:digit:]]+)*");
  regex intervalPat("(\\-?[[:digit:]]+)\\-(\\-?[[:digit:]]+)");
  std::smatch sm;

  const string resiInput = input.substr(5);

  if( regex_match(resiInput, intListPat) ){
    vector<string> tokens = Util::split(resiInput, '+');
    for(const string& token: tokens){
      int resi = std::stoi(token);
      m_residueIDs.insert(resi);
    }
  }else if(regex_match(resiInput, sm, intervalPat)){
    int intervalStart = stoi(sm[1]);
    int intervalEnd = stoi(sm[2]);
    for(int i=std::min(intervalStart, intervalEnd); i<=std::max(intervalStart, intervalEnd); i++){
      m_residueIDs.insert(i);
    }
  }else{
    throw std::runtime_error("KGS error: Selection::ResiClause - Can't parse \""+input+"\"");
  }
}
bool Selection::ResiClause::inSelection(const Atom* a) const{
  return m_residueIDs.find(a->getResidue()->getId()) != m_residueIDs.end();
}

Selection::ResnClause::ResnClause(const std::string& input)
{
  regex strListPat("([^\\+ ])+(\\+[^+ ]+)*");
  const string resnInput = input.substr(5);

  if( regex_match(resnInput, strListPat) ){
    vector<string> tokens = Util::split(resnInput, '+');
    for(const string& token: tokens){
      m_residueNames.push_back(token);
    }
  }
}
bool Selection::ResnClause::inSelection(const Atom* a) const{
  return std::find( m_residueNames.begin(), m_residueNames.end(), a->getResidue()->getName() )
         != m_residueNames.end();
}

Selection::NameClause::NameClause(const std::string& input)
{
  regex strListPat("([^\\+ ])+(\\+[^+ ]+)*");
  const string nameInput = input.substr(5);

  if( regex_match(nameInput, strListPat) ){
    vector<string> tokens = Util::split(nameInput, '+');
    for(const string& token: tokens){
      m_atomNames.push_back(token);
    }
  }
}
bool Selection::NameClause::inSelection(const Atom* a) const
{
  return std::find( m_atomNames.begin(), m_atomNames.end(), a->getName() )
         != m_atomNames.end();
}

Selection::ElemClause::ElemClause(const std::string& input)
{
  regex strListPat("([^\\+ ])+(\\+[^+ ]+)*");
  const string elemInput = input.substr(5);

  if( regex_match(elemInput, strListPat) ){
    vector<string> tokens = Util::split(elemInput, '+');
    for(const string& token: tokens){
      m_atomElements.push_back(token);
    }
  }
}
bool Selection::ElemClause::inSelection(const Atom* a) const
{
  return std::find( m_atomElements.begin(), m_atomElements.end(), a->getElement() )
         != m_atomElements.end();
}

Selection::AllClause::AllClause(const std::string& input){}
bool Selection::AllClause::inSelection(const Atom* a) const{
  return true;
}

Selection::BackboneClause::BackboneClause(const std::string& input):
    NameClause("name N+CA+C+P+O5'+C5'+C4'+C3'+O3'")
{}

Selection::HeavyClause::HeavyClause(const std::string& input):
    NotClause("not hydro")
{}

Selection::HydroClause::HydroClause(const std::string& input):
    ElemClause("elem H")
{}

std::ostream& operator<<(std::ostream& os, const Selection& s)
{
  return os<<s.m_selectionPattern;
}

std::ostream& operator<<(std::ostream& os, const Selection* s)
{
  return os<<s->m_selectionPattern;
}

//Selection::Selection( ) : delim_( " " ) {
//}
////---------------------------------------------------------
//Selection::Selection( string selection ) : delim_( " " ),
//					   selection_( selection ) {
//	split( selection_, delim_, selectionWords_ );
//}
////---------------------------------------------------------
//Selection::Selection( string selection, string delim ) : delim_( delim ),
//					       		 selection_( selection ) {
//	split( selection_, delim_, selectionWords_ );
//}
////---------------------------------------------------------
//void Selection::print() const {
//	cout<<"Selection)\t selection = \""<<selection()<<"\""<<endl;
//	cout<<"Selection)\t selectionWords = ";
//	for(vector<string>::const_iterator it = selectionWords_.begin(); it != selectionWords_.end(); ++it) {
//		cout<<"\'"<<(*it);
//		//if( it != prev(selectionWords_.end()) ) cout<<"\', "; //c++11
//		if( (it+1) != selectionWords_.end() ) cout<<"\', ";
//		else cout<<"\'";
//	}
//	cout<<endl;
//}
////---------------------------------------------------------
//void Selection::print( string selName ) const {
//	cout<<"Selection)\t"<<selName<<"  selection = \""<<selection()<<"\""<<endl;
//	cout<<"Selection)\t"<<selName<<"  selectionWords = ";
//	for(vector<string>::const_iterator it = selectionWords_.begin(); it != selectionWords_.end(); ++it) {
//		cout<<"\'"<<(*it);
//		//if( it != prev(selectionWords_.end()) ) cout<<"\', "; //c++11
//		if( (it+1) != selectionWords_.end() ) cout<<"\', ";
//		else cout<<"\'";
//	}
//	cout<<endl;
//}
////---------------------------------------------------------
//// Accessor
//string Selection::delim() const {
//	return delim_;
//}
////---------------------------------------------------------
//// Mutator
//void Selection::delim(string delim) {
//	delim_ = delim;
//}
//
////---------------------------------------------------------
//// Accessor
//string Selection::selection() const {
//	return selection_;
//}
////---------------------------------------------------------
//// Mutator
//void Selection::selection(string selection) {
//	selection_ = selection;
//	// The selectionWords_ vector of strings has to contain the new content
//	split( selection_, delim_, selectionWords_ );
//}
////---------------------------------------------------------
//// Accessor
//vector<string> Selection::selectionWords( ) const {
//	return selectionWords_;
//}
////---------------------------------------------------------
//// Mutator
//void Selection::selectionWords( vector<string> selectionWords ) {
//	selectionWords_ = selectionWords;
//	// The selection_ string has to be made of the new words
//	combine( selectionWords_, delim_, selection_ );
//}
////---------------------------------------------------------
///*
// * Parses a [RESIDUE ID] selection and returns a vector of corresponding residues in m_molecule.
// * Assumes the selection starts with either "resid" or "not" followed by [ranges of] resids ("resid 1 to 123 145 200 300 to 400").
// * TODO: Rough parsing, need to develop a more elegant routine (keywords, boolean and/or/not, etc).
// */
//vector<Residue*>
//Selection::getSelectedResidues( const Molecule *protein ) const {
//	vector<Residue*> residues;
//	vector<string> words = selectionWords();
//	bool isInSelection;
//
//	if( (*words.begin()) == "resid" ||
//			//( (*words.begin()) == "not" && *next(words.begin()) == "resid" ) ) {
//			( (*words.begin()) == "not" && *(words.begin()+1) == "resid" ) ) {
//		vector<Residue*> residuesProtein = protein->chains[0]->getResidues();
//		unsigned int res1ID, res2ID = 0;
//		//vector<string>::const_iterator itwstart = next(words.begin());
//		vector<string>::const_iterator itwstart = (words.begin()+1);
//		if( (*words.begin()) == "not" )
//			itwstart = (words.begin()+2);
//			//itwstart = next(words.begin(),2);
//
//		for ( vector<Residue*>::const_iterator itres = residuesProtein.begin();
//						    	  itres != residuesProtein.end(); ++itres ) {
//			Residue* res = *itres;
//			isInSelection = false;
//			res2ID=0;
//			for(vector<string>::const_iterator itw = itwstart; itw != words.end(); ++itw)
//				if( (*itw) != "to" ) {
//					res1ID = atoi( (*itw).c_str() );
//					if ( res2ID == res1ID ) continue;
//					//if ( next(itw) != words.end() ) {
//					if ( (itw+1) != words.end() ) {
//						if ( (*(itw+1)) != "to" || (itw+1) == words.end() )
//						//if ( (*next(itw)) != "to" || next(itw) == words.end() )
//							res2ID = res1ID;
//						else
//							res2ID = atoi ( ( *(itw+2) ).c_str() );
//							//res2ID = atoi ( ( *(next(itw,2) ).c_str() );
//					}else
//						res2ID = res1ID;
//
//					if( res->isWithinResidueRange( res1ID, res2ID ) ) {
//						if( (*words.begin()) != "not" ) residues.push_back( res );
//						else isInSelection = true;
//						break;
//					}
//				}
//			if( (*words.begin()) == "not" && !isInSelection ) residues.push_back( res );
//		}
//	}
//	return residues;
//}
////---------------------------------------------------------
///*
// * Parses a [ATOM NAME] selection and returns a vector of corresponding atoms in m_molecule.
// * Assumes the selection starts with either "name" or "not" followed by atom names ("name CA N").
// * TODO: Rough parsing, need to develop a more elegant routine (keywords, boolean and/or/not, etc).
// */
//vector<Atom*>
//Selection::getSelectedAtoms( const Molecule *protein )
//{
//	vector<Residue*> residues = protein->chains[0]->getResidues();
//	return getSelectedAtoms( residues );
//}
////---------------------------------------------------------
///*
// * Parses a [ATOM NAME] selection and returns a vector of corresponding atoms in input m_molecule residues.
// * Assumes the selection starts with either "name" or "not" followed by atom names ("name CA N").
// * TODO: Rough parsing, need to develop a more elegant routine (keywords, boolean and/or/not, etc).
// */
//vector<Atom*>
//Selection::getSelectedAtoms( const vector<Residue*> residues ) {
//	vector<Atom*> atoms;
//	vector<string> words = selectionWords();
//	bool isInSelection;
//
//	// TODO: This has to go to a separate routine/function
//	if( (*words.begin()) == "backbone" ) selection( "name N CA C O" );
//	else if( (*words.begin()) == "not" && *(words.begin()+1) == "backbone" ) selection( "not name N CA C O" );
//	//else if( (*words.begin()) == "not" && *next(words.begin()) == "backbone" ) selection( "not name N CA C O" );
//	else if( (*words.begin()) == "heavy" ) selection( "not name H*" );
//	else if( (*words.begin()) == "not" && *(words.begin()+1) == "heavy" ) selection( "name H*" );
//	//else if( (*words.begin()) == "not" && *next(words.begin()) == "heavy" ) selection( "name H*" );
//	else if( (*words.begin()) == "all" ) selection( "name C* N* O* P* S* H*" );
//	else if( (*words.begin()) == "none" ) return atoms;
//
//	words = selectionWords();
//
//	if( (*words.begin()) == "name" ||
//			( (*words.begin()) == "not" && *(words.begin()+1) == "name" ) ) {
//			//( (*words.begin()) == "not" && *next(words.begin()) == "name" ) ) {
//		vector<string>::const_iterator itwstart = (words.begin()+1);
//		//vector<string>::const_iterator itwstart = next(words.begin());
//		if( (*words.begin()) == "not" )
//			itwstart = (words.begin()+2);
//			//itwstart = next(next(words.begin()));
//		for ( vector<Residue*>::const_iterator itres = residues.begin();
//						    	  itres != residues.end(); ++itres ) {
//			Residue* res = *itres;
//			if( (*words.begin()) != "not" ) {
////				//Todo: this does not work for names with * in the end --> case "all", case "not heavy"
////				for(vector<string>::const_iterator itw = itwstart; itw != words.end(); ++itw)
////					if( res->Atom_map_by_name.find(*itw) != res->Atom_map_by_name.end() ){
////						cout<<"Pushing back: "<<res->Atom_map_by_name.find(*itw)->second<<endl;
////						atoms.push_back( res->Atom_map_by_name.find(*itw)->second );
////					}
//				// Loop over all the atoms of the residue
//				for(auto const& atom: res->getAtoms() ) {
//					isInSelection = false;
//					for(vector<string>::const_iterator itw = itwstart; itw != words.end(); ++itw)
//						if( (*itw) == "H*" || (*itw) == "C*" || (*itw) == "N*" ||
//						    		      (*itw) == "O*" || (*itw) == "P*" || (*itw) == "S*" ) {
//							if( atom->compareType( (*itw).substr(0,1) ) ) {
//								isInSelection = true;
//								break;
//							}
//						}else{
//							if( atom->compareName(*itw) ) {
//								isInSelection = true;
//								break;
//							}
//						}
//          if( isInSelection ) atoms.push_back( atom );
//				}
//			}else{
//				// Loop over all the atoms of the residue
//        for(auto const& atom: res->getAtoms() ) {
//					isInSelection = false;
//					for(vector<string>::const_iterator itw = itwstart; itw != words.end(); ++itw)
//						if( (*itw) == "H*" || (*itw) == "C*" || (*itw) == "N*" ||
//						    		      (*itw) == "O*" || (*itw) == "P*" || (*itw) == "S*" ) {
//              if( atom->compareType( (*itw).substr(0,1) ) ) {
//								isInSelection = true;
//								break;
//							}
//						}else{
//              if( atom->compareName(*itw) ) {
//								isInSelection = true;
//								break;
//							}
//						}
//					if( !isInSelection ) atoms.push_back( atom );
//				}
//			}
//		}
//	}
//	return atoms;
//}
//
////---------------------------------------------------------
///*
// * Combine a vector of strings.
// * Result returned in a pre-constructed string, separated by delimiter.
// */
//string& Selection::combine( const vector<string> &words,
//		       	    string delim,
//		            string &s ) {
//	s.clear();
//	for(vector<string>::const_iterator it = words.begin(); it != words.end(); ++it) {
//		s.append( (*it) );
//		if( it != (words.end()-1) ) s.append( delim_ );
//		//if( it != prev(words.end()) ) s.append( delim_ ); //c++11
//	}
//	return s;
//}
////---------------------------------------------------------
///*
// * Combine a vector of strings.
// * Result returned in a new string, separated by delimiter.
// */
//string Selection::combine( const vector<string> &words,
//		           string delim ) {
//	string s;
//	combine(words, delim, s);
//	return s;
//}
////---------------------------------------------------------
///*
// * Split string by delimiter (uses 1st char).
// * Results returned in a pre-constructed vector.
// */
//vector<string>& Selection::split( const string &s,
//		       	          string delim,
//		       		  vector<string> &words ) {
//	words.clear();
//	stringstream ss(s);
//	string item;
//	char c = delim[0]; // assumes 1st char of delim string is the separator!
//	while( getline(ss, item, c) )
//	        words.push_back(item);
//	return words;
//}
////---------------------------------------------------------
///*
// * Split string by delimiter (uses 1st char).
// * Results returned in a new vector.
// */
//vector<string> Selection::split( const string &s,
//		      		 string delim ) {
//	vector<string> words;
//	split(s, delim, words);
//	return words;
//}
////---------------------------------------------------------
///*
// * Split string by delimiter (uses 1st char).
// * Results returned in a pre-constructed vector of ints.
// */
//vector<int>& Selection::split( const string &s,
//		       	          string delim,
//		       		  vector<int> &numbers ) {
//	numbers.clear();
//	stringstream ss(s);
//	string item;
//	char c = delim[0]; // assumes 1st char of delim string is the separator!
//	while( getline(ss, item, c) ){
//		int itemNum = atoi(item.c_str());
//		numbers.push_back(itemNum);
//	}
//	return numbers;
//}


