
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "core/Molecule.h"
#include "core/Chain.h"
#include "core/Grid.h"
#include "core/Atom.h"
#include "JacobianRelated.h"
#include "Util.h"
#include "HbondIdentifier.h"
#include "IO.h"
#include "math/MathUtility.h"
#include "DisjointSets.h"
#include "core/ProteinHBond.h"
#include "Logger.h"


using namespace std;

/*
bool sugar_valid2(Vector3& p0, Vector3& p1, Vector3& p2, Vector3& p3, Vector3& p4, double Am){
	double torsion = Am;

	double angles[5];
	double dists[5];
	angles[0] = VectorAngle(p1-p0,p4-p0);
	angles[1] = VectorAngle(p2-p1,p0-p1);
	angles[2] = VectorAngle(p3-p2,p1-p2);
	angles[3] = VectorAngle(p4-p3,p2-p3);
	angles[4] = VectorAngle(p0-p4,p3-p4);
	dists[0] = (p0-p1).length();
	dists[1] = (p1-p2).length();
	dists[2] = (p2-p3).length();
	dists[3] = (p3-p4).length();
	dists[4] = (p4-p0).length();

	Vector3 p0x = p1-p0; normalize(p0x);
	Vector3 p0z = cross(p0x, p4-p0); normalize(p0z);
	Vector3 p0y = cross(p0z, p0x); 
	Vector3 p2_ = p1+(cos(CTK_PI-angles[1])*dists[1])*p0x + (sin(CTK_PI-angles[1])*cos(torsion)*dists[1])*p0y + (sin(CTK_PI-angles[1])*sin(torsion)*dists[1])*p0z;

	Vector3 n = (p2_-p1)/dists[1];
	Vector3 S1 = p2_+n*( cos(CTK_PI-angles[2])*dists[2] );
	double S2off = n.dot(S1-p4);
	Vector3 S2 = p4 + n*S2off;
	double l = (S1-S2).length();
	double r1 = sin(CTK_PI-angles[2])*dists[2];
	double r2 = ( dists[3]*dists[3] - S2off*S2off );//Still need to take the square m_root
	if( r2<0 ) return false;
	r2 = sqrt(r2);
	if(r1+r2<l) return false;
	return true;
}
 */


int main( int argc, char* argv[] ) {
  enableLogger("torsions");

  log("torsions") << setw(11) << "#molecule" << " ";
  log("torsions") << setw(4) << "res" << " ";
  log("torsions") << setw(9) << "alpha" << " ";
  log("torsions") << setw(9) << "beta" << " ";
  log("torsions") << setw(9) << "gamma" << " ";
  log("torsions") << setw(9) << "delta" << " ";
  log("torsions") << setw(9) << "epsilon" << " ";
  log("torsions") << setw(9) << "zeta" << " ";
  log("torsions") << setw(9) << "tau0" << " ";
  log("torsions") << setw(9) << "tau1" << " ";
  log("torsions") << setw(9) << "tau2" << " ";
  log("torsions") << setw(9) << "tau3" << " ";
  log("torsions") << setw(9) << "tau4" << " ";
  log("torsions") << setw(9) << "chi" << " ";
  log("torsions") << setw(9) << "up" << " ";
  log("torsions") << setw(9) << "Am" << " ";
  log("torsions") << setw(9) << "tau" << " ";
  log("torsions") << endl;

  for (int a = 1; a < argc; a++) {
    char *tmp = realpath(argv[a], nullptr);
    if (tmp == nullptr) {
      cerr << argv[a] << " is not a valid PDB-file" << endl;
      exit(-1);
    }
    string pdb_file(tmp);
    int nameSplit = pdb_file.find_last_of("/\\");
    string protein_name = pdb_file.substr(nameSplit + 1);
    int pos = protein_name.rfind(".pdb");
    if (pos != string::npos) protein_name = protein_name.substr(0, pos);

    vector<string> extraCovBonds;
    Molecule *protein = new Molecule();
    IO::readPdb(protein, argv[a], extraCovBonds);

    IO::readRigidbody(protein);
    protein->buildSpanningTree();

    bool rna = true;
    for (int i = 0; i < protein->atoms.size(); i++) {
      if (protein->atoms[i]->getName() == "CA") {
        rna = false;
        break;
      }
    }

    for (auto const &chainPtr: protein->chains) {
      string chain = chainPtr->getName();
      for (int r = protein->getMinResidueNumber(); r <= protein->getMaxResidueNumber(); r++) {
        if (rna) {
          Atom *O3prev = protein->getAtom(chain, r - 1, "O3'");
          Atom *P = protein->getAtom(chain, r, "P");
          Atom *O5 = protein->getAtom(chain, r, "O5'");
          Atom *C5 = protein->getAtom(chain, r, "C5'");
          Atom *C4 = protein->getAtom(chain, r, "C4'");
          Atom *C3 = protein->getAtom(chain, r, "C3'");
          Atom *C2 = protein->getAtom(chain, r, "C2'");
          Atom *C1 = protein->getAtom(chain, r, "C1'");
          Atom *O4 = protein->getAtom(chain, r, "O4'");
          Atom *O3 = protein->getAtom(chain, r, "O3'");
          Atom *Pnext = protein->getAtom(chain, r + 1, "P");
          Atom *O5next = protein->getAtom(chain, r + 1, "O5'");

          Atom *scN = protein->getAtom(chain, r, "N9");
          Atom *scC = nullptr;
          if (scN != nullptr) {//purine
            scC = protein->getAtom(chain, r, "C4");
          } else {//pyrimidine
            scN = protein->getAtom(chain, r, "N1");
            scC = protein->getAtom(chain, r, "C2");
          }

          double alpha = (O3prev == nullptr || P == nullptr) ? -42.0 : TorsionalAngle(O3prev->m_Position, P->m_Position,
                                                                                O5->m_Position, C5->m_Position);
          double beta = (P == nullptr) ? -42.0 : TorsionalAngle(P->m_Position, O5->m_Position, C5->m_Position, C4->m_Position);
          double gamma = TorsionalAngle(O5->m_Position, C5->m_Position, C4->m_Position, C3->m_Position);
          double delta = TorsionalAngle(C5->m_Position, C4->m_Position, C3->m_Position, O3->m_Position);
          double epsilon = (Pnext == nullptr) ? -42.0 : TorsionalAngle(C4->m_Position, C3->m_Position, O3->m_Position,
                                                                    Pnext->m_Position);
          double zeta = (Pnext == nullptr || O5next == nullptr) ? -42.0 : TorsionalAngle(C3->m_Position, O3->m_Position,
                                                                                   Pnext->m_Position, O5next->m_Position);

          double tau0 = TorsionalAngle(C4->m_Position, O4->m_Position, C1->m_Position, C2->m_Position);
          double tau1 = TorsionalAngle(O4->m_Position, C1->m_Position, C2->m_Position, C3->m_Position);
          double tau2 = TorsionalAngle(C1->m_Position, C2->m_Position, C3->m_Position, C4->m_Position);
          double tau3 = TorsionalAngle(C2->m_Position, C3->m_Position, C4->m_Position, O4->m_Position);
          double tau4 = TorsionalAngle(C3->m_Position, C4->m_Position, O4->m_Position, C1->m_Position);
          double chi = TorsionalAngle(O4->m_Position, C1->m_Position, scN->m_Position, scC->m_Position);

          log("torsions") << setw(11) << protein_name << " ";
          log("torsions") << setw(4) << C1->getResidue()->getId() << " ";
          log("torsions") << setw(9) << setprecision(4) << alpha << " ";
          log("torsions") << setw(9) << setprecision(4) << beta << " ";
          log("torsions") << setw(9) << setprecision(4) << gamma << " ";
          log("torsions") << setw(9) << setprecision(4) << delta << " ";
          log("torsions") << setw(9) << setprecision(4) << epsilon << " ";
          log("torsions") << setw(9) << setprecision(4) << zeta << " ";
          log("torsions") << setw(9) << setprecision(4) << tau0 << " ";
          log("torsions") << setw(9) << setprecision(4) << tau1 << " ";
          log("torsions") << setw(9) << setprecision(4) << tau2 << " ";
          log("torsions") << setw(9) << setprecision(4) << tau3 << " ";
          log("torsions") << setw(9) << setprecision(4) << tau4 << " ";
          log("torsions") << setw(9) << setprecision(4) << chi << " ";
          //log("torsions")<<setw(9)<<setprecision(4)<<vs->initUp;
          //log("torsions")<<setw(9)<<setprecision(4)<<vs->Amplitude;
          //log("torsions")<<setw(9)<<setprecision(4)<<acos(vs->initTorsion/vs->Amplitude)*vs->initUp;
          log("torsions") << endl;
        } else {
          Atom *Cprev = protein->getAtom(chain, r - 1, "C");
          Atom *N = protein->getAtom(chain, r, "N");
          Atom *CA = protein->getAtom(chain, r, "CA");
          Atom *C = protein->getAtom(chain, r, "C");
          Atom *Nnext = protein->getAtom(chain, r + 1, "N");
          Atom *CAnext = protein->getAtom(chain, r + 1, "CA");

          double phi = (Cprev == nullptr || N == nullptr) ? 180 : TorsionalAngle(Cprev->m_Position, N->m_Position, CA->m_Position,
                                                                           C->m_Position);
          double psi = (N == nullptr || Nnext == nullptr) ? 180 : TorsionalAngle(N->m_Position, CA->m_Position, C->m_Position,
                                                                           Nnext->m_Position);
          double omega = (CA == nullptr || CAnext == nullptr) ? 180 : TorsionalAngle(CA->m_Position, C->m_Position,
                                                                               Nnext->m_Position, CAnext->m_Position);

          log("torsions") << setw(11) << protein_name << " ";
          log("torsions") << setw(4) << CA->getResidue()->getId() << " ";
          log("torsions") << setw(9) << setprecision(4) << phi << " ";
          log("torsions") << setw(9) << setprecision(4) << psi << " ";
          log("torsions") << setw(9) << setprecision(4) << omega << " ";
          log("torsions") << endl;

        }
      }
    }
  }
}


