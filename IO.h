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
#ifndef IO_H
#define IO_H

#include <string>
#include <vector>

#include "core/Molecule.h"
#include "core/Residue.h"
#include "Util.h"
#include "Selection.h"
#include "ResidueProfiles.h"


class IO {
  public:
	static ResidueProfile readResidueProfile ();
	static void readPdb (Molecule * protein, std::string pdb_file, std::vector<std::string> &hbondsAsCov, Molecule * reference = nullptr);
	static void readDssp (Molecule * protein, std::string dssp_file);
	static void readRigidbody (Molecule * molecule);
  static void readRigidbody (Molecule * molecule, Selection& movingResidues);
	static void writePdb (Molecule * molecule, std::string output_file_name);
	static void writePyMolScript(Molecule * molecule, std::string pdb_file, std::string output_file_name);
  static void writeBondLengthsAndAngles (Molecule *molecule, std::string output_file_name);
	static void writeCovBonds (Molecule *molecule, std::string output_file_name);
	static void readCovBonds  (Molecule *molecule, std::string input_file_name);
	static void writeHbonds (Molecule * molecule, std::string output_file_name);
	static void writeHbondsChange (Molecule * molecule, std::string output_file_name);
	static void readHbonds (Molecule *molecule, std::string hbond_file_name);
  static void readAnnotations (Molecule *molecule, std::string annotation_file_name);
	static void readHbonds_dssr(Molecule * molecule, std::string dssrFile);
	static void readHbonds_rnaview(Molecule * molecule, std::string file, bool fillAnnotations);
	static void readHbonds_first(Molecule * molecule, std::string file);
	static void readHbonds_vadar(Molecule * molecule, std::string file);
	static void writeRBs(Molecule * molecule, std::string output_file_name);
	static void writeStats(Molecule * molecule, std::string output_file_name);
	static void writeQ (Molecule *molecule, Configuration* referenceConf, std::string output_file_name);
  static void writeTrajectory (Molecule *molecule, std::string output_file, std::string output_mdl, Molecule *target = nullptr);
  private:
	static void makeCovBond (Residue* res1, Residue* res2, std::string atom_name1, std::string atom_name2);
};

#endif
