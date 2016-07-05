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
#ifndef SELECTION_H
#define SELECTION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Residue.h"


/**
 * A selection of atoms in a molecule based on a pymol-like selection-pattern. The
 * selection-pattern has the following grammar, where the top-most operators bind the
 * weakest:
 * <pre>
 * SEL :== SEL + SEL         # Disjunction of two selections
 *     :== SEL or SEL        # Disjunction of two selections
 *     :== SEL & SEL         # Conjunction of two selections
 *     :== SEL and SEL       # Conjunction of two selections
 *     :== not SEL           # Negation of selection
 *     :== resi <int>-<int>  # Residue ID interval
 *     :== resi INTLIST      # Residue ID list
 *     :== resn STRINGLIST   # Residue name list
 *     :== name STRINGLIST   # Atom name list
 *     :== elem STRINGLIST   # Element name list
 *     :== all               # All atoms selected
 *     :== backbone          # Backbone atoms only: "name CA+C+N+P+O5'+C5'+C4'+C3'+O3'"
 *     :== heavy             # Heavy atoms only: "not elem H"
 *     :== hydro             # Hydrogen atoms only: "elem H"
 * INTLIST    :== <int>
 *            :== <int>+INTLIST
 * STRINGLIST :== <string>
 *            :== <string>+STRINGLIST
 * </pre>
 *
 * There are two uses of the symbol '+' and the whitespace around it is important. A
 * typical selection will be composed of a series of clauses or-ed together, e.g.
 * "resi 5 or resi 10" will select all atoms in residues with id 5 and 10.
 *
 * Since "and" binds stronger than "or", the pattern "name CA and resi 2 or resi 3" will
 * select only the C-alpha atom of residue 2 but also all atoms in residue 3.
 *
 * There is no parenthesis in selection patterns, so to select all C-alphas in two loop
 * regions use the selection "name CA and resi 10-15 + name CA and resi 20-25".
 *
 * To use a selection object, first instantiate it with the pattern and then call either
 * getSelectedAtoms, getSelectedBonds, or getSelectedResidues on a molecule object. The
 * result will be cached so subsequent calls on the same molecule object will be fast.
 * <pre>
 * Selection sel("resi 10-15 and backbone");
 * Molecule mol;
 * //Initialize mol
 * //...
 * for( auto atom: sel.getSelectedAtoms(&mol))
 *   cout<<atom<<endl;
 * </pre>
 *
 * Passing an empty string corresponds to the pattern "all".
 */
class Selection {
 public:

  /** Equivalent to constructing with the string "all" */
  Selection();

  /** Initialize the selection with the specified selection-pattern. */
  Selection( const std::string& pattern );

  /** Return all atoms in mol that match the selection-pattern passed to the constructor. The first time
   * this function is called it will iterate through all atoms in mol and the result will be cached.
   * Subsequent calls will simply return a reference to the vector computed in the first call. */
  std::vector<Atom *>& getSelectedAtoms( const Molecule *mol );

  /** Return all bonds for which end-point atoms are matched by the selection pattern. The first time
   * this function is called it will iterate through all bonds in mol and the result will be cached.
   * Subsequent calls will simply return a reference to the vector computed in the first call.  */
  std::vector<Bond *>& getSelectedBonds( const Molecule *mol );

  /** Return all residues for which all contained atoms are matched by the selection pattern. The first time
   * this function is called it will iterate through all atoms in mol and the result will be cached.
   * Subsequent calls will simply return a reference to the vector computed in the first call.  */
  std::vector<Residue*>& getSelectedResidues( const Molecule *mol );


 private:
  class Clause; //Forward decl

  const std::string m_selectionPattern;

  /// The parsed representation of the selectionPattern
  const Clause* m_rootClause;

  std::map<const Molecule *, std::vector<Atom *> > m_cachedAtoms;
  std::map<const Molecule *, std::vector<Bond *> > m_cachedBonds;
  std::map<const Molecule *, std::vector<Residue *> > m_cachedResidues;

  static Clause* parseClause(const std::string& input);

  class Clause{
   public:
    virtual bool inSelection(Atom* a) const = 0;
  };

  class OrClause: public Clause{
   public:
    OrClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::vector<Clause*> m_childClauses;
  };

  class AndClause: public Clause{
   public:
    AndClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::vector<Clause*> m_childClauses;
  };

  class NotClause: public Clause{
   public:
    NotClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    const Clause* m_childClause;
  };

  class ResiClause: public Clause{
   public:
    ResiClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::set<int> m_residueIDs;
  };

  class ResnClause: public Clause{
   public:
    ResnClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::vector<std::string> m_residueNames;
  };

  class NameClause: public Clause{
   public:
    NameClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::vector<std::string> m_atomNames;
  };

  class ElemClause: public Clause{
   public:
    ElemClause(const std::string& input);
    bool inSelection(Atom* a) const;
   private:
    std::vector<std::string> m_atomElements;
  };

  class AllClause: public Clause{
   public:
    AllClause(const std::string& input);
    bool inSelection(Atom* a) const;
  };

  class BackboneClause: public NameClause{
   public:
    BackboneClause(const std::string& input);
  };

  class HeavyClause: public NotClause{
   public:
    HeavyClause(const std::string& input);
  };

  class HydroClause: public ElemClause{
   public:
    HydroClause(const std::string& input);
  };

};

//public:
//  Selection( );
//	Selection( std::string selection );
//	Selection( std::string selection, std::string delim );
//	void print() const;
//	void print( std::string selName ) const;
//
//	// Mutator and Accessor
//	void delim(std::string delim);
//	std::string delim() const;
//	// Mutator and Accessor
//	void selection(std::string selection);
//	std::string selection() const;
//	// Mutator and Accessor
//	void selectionWords( std::vector<std::string> selectionWords );
//	std::vector<std::string> selectionWords() const;
//
//	std::vector<Residue*> getSelectedResidues( const Molecule *protein ) const;
//	std::vector<Atom*> getSelectedAtoms( const Molecule *protein );
//	std::vector<Atom*> getSelectedAtoms( const std::vector<Residue*> residues );
//
//	static std::vector<std::string> &split( const std::string &s, std::string delim, std::vector<std::string> &words );
//	static std::vector<std::string>  split( const std::string &s, std::string delim );
//	static std::vector<int> &split( const std::string &s, std::string delim, std::vector<int> &numbers );
//	std::string &combine( const std::vector<std::string> &words, std::string delim, std::string &s );
//	std::string  combine( const std::vector<std::string> &words, std::string delim );
//
//private:
//	std::string selection_, delim_;
//	std::vector<std::string> selectionWords_;
//};

#endif // SELECTION_H

