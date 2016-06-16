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
#include "ProteinHBond.h"
#include "core/Atom.h"
#include "math/MathUtility.h"
#include "Logger.h"

using namespace Math3D;
using namespace std;

Hbond::Hbond(Atom* hatom, Atom* acceptor, Atom* donor, Atom* aa, double energy) {
	Hatom = hatom;
	Acceptor = acceptor;
	Donor = donor;
	AA = aa;

	// super-class attributes
	Atom1 = hatom;
	Atom2 = acceptor;
	BondType = "HB";
	Bars = 5;
	constrained = false;
	m_iniDist_H_A = VectorLength(Atom1->m_Position,Atom2->m_Position);
	m_iniAngle_D_H_A = getAngle_D_H_A();
	m_iniAngle_H_A_AA = getAngle_H_A_AA();

  identifyHybridization();

  m_iniEnergy = energy;
  if(energy == DEFAULT_HBOND_ENERGY)
    m_iniEnergy = computeEnergy();

//  Vector3 x,y,z;
//  coordinateSystem(hatom, x,y,z);
//  Vector3 d = acceptor->m_Position -hatom->m_Position;
//  idealA[0] = x.dot(d);
//  idealA[1] = y.dot(d);
//  idealA[2] = z.dot(d);
//  coordinateSystem(acceptor, x,y,z);
//  idealH[0] = -x.dot(d);
//  idealH[1] = -y.dot(d);
//  idealH[2] = -z.dot(d);
}

Hbond::Hbond(Hbond & hbond) {
	Hatom = hbond.Hatom;
	Acceptor = hbond.Acceptor;
	Donor = hbond.Donor;
	AA = hbond.AA;
	m_iniEnergy = hbond.m_iniEnergy;

	// super-class attributes
	Atom1 = hbond.Atom1;
	Atom2 = hbond.Atom2;
	BondType = hbond.BondType;
	Bars = hbond.Bars;
	constrained = hbond.constrained;
	m_iniDist_H_A = hbond.m_iniDist_H_A;
	m_iniAngle_D_H_A = hbond.m_iniAngle_D_H_A;
	m_iniAngle_H_A_AA = hbond.m_iniAngle_H_A_AA;

  m_D_sp2 = hbond.m_D_sp2;
  m_D_sp3 = hbond.m_D_sp3;
  m_A_sp2 = hbond.m_A_sp2;
  m_A_sp3 = hbond.m_A_sp3;

//	idealA = hbond.idealA;
//	idealH = hbond.idealH;
}


bool Hbond::isSame (Hbond * b2) {
	if ( Hatom->getName() == b2->Hatom->getName() &&
		Hatom->getResidue()->getId() == b2->Hatom->getResidue()->getId() &&
		Donor->getName()==b2->Donor->getName() &&
		Donor->getResidue()->getId() == b2->Donor->getResidue()->getId() &&
		Acceptor->getName()==b2->Acceptor->getName() &&
		Acceptor->getResidue()->getId() == b2->Acceptor->getResidue()->getId() &&
		AA->getName()==b2->AA->getResidue()->getName() &&
		AA->getResidue()->getId() == b2->AA->getResidue()->getId() ){
			return true;
	}
	return false;
}

//Vector3 Hbond::getIdealHPoint(){
//    Vector3 x,y,z;
//    coordinateSystem(Acceptor, x,y,z);
//    return Acceptor->m_Position + (x*idealH[0]) + (y*idealH[1]) + (z*idealH[2]);
//}
//
//Vector3 Hbond::getIdealAcceptorPoint(){
//    Vector3 x,y,z;
//    coordinateSystem(Hatom, x,y,z);
//    return Hatom->m_Position + (x*idealA[0]) + (y*idealA[1]) + (z*idealA[2]);
//}

//void Hbond::coordinateSystem(Atom* a, Vector3& x, Vector3& y, Vector3& z ){
//	Atom* a1 = a;
//    Atom* a2 = a->Cov_neighbor_list[0];
//    Atom* a3 = nullptr;
//    if( a1->Cov_neighbor_list.size()>1) a3 = a1->Cov_neighbor_list[1];
//    else{
//        if(a2->Cov_neighbor_list[0]!=a1)
//            a3 = a2->Cov_neighbor_list[0];
//        else
//            a3 = a2->Cov_neighbor_list[1];
//    }
//
//    x = a1->m_Position -a2->m_Position;
//    y = a3->m_Position -a2->m_Position;
//    z = cross(x,y);
//}

double Hbond::getLength() {
  double length = Hatom->m_Position.distanceTo(Acceptor->m_Position);
  return length;
}

double Hbond::getDistance_D_A() {
  double distance = Donor->m_Position.distanceTo(Acceptor->m_Position);
  return distance;
}

double Hbond::getAngle_D_H_A() {
  double ret = Angle(Donor->m_Position, Hatom->m_Position, Acceptor->m_Position);
  return ret;
}

double Hbond::getAngle_H_A_AA() {
  double ret = Angle(Hatom->m_Position, Acceptor->m_Position, AA->m_Position);
  return ret;
}

double Hbond::getOutOfPlaneAngle() {

	//Todo: Check if one should rather use the other covalent neighbors instead of H-atom and Acceptor

	Atom *a1, *a2, *a3;

  if(Donor->Cov_neighbor_list.size() >= 3) {
    a2 = Donor->Cov_neighbor_list[0];
    a3 = Donor->Cov_neighbor_list[1];
    if (a2 == Hatom)
      a2 = Donor->Cov_neighbor_list[2];
    if (a3 == Hatom)
      a3 = Donor->Cov_neighbor_list[2];
    a1 = Hatom;
  }
  else{
    a1=Donor;
    a2 = Hatom;
    a3 =  Donor->Cov_neighbor_list[0] == Hatom ? Donor->Cov_neighbor_list[1] : Donor->Cov_neighbor_list[0];
  }

  Math3D::Vector3 normal1 = UnitNormal(a1->m_Position,a2->m_Position, a3->m_Position);

  if(AA->Cov_neighbor_list.size() >= 3) {
    a2 = AA->Cov_neighbor_list[0];
    a3 = AA->Cov_neighbor_list[1];
    if (a2 == Acceptor)
      a2 = AA->Cov_neighbor_list[2];
    if (a3 == Acceptor)
      a3 = AA->Cov_neighbor_list[2];
    a1 = Acceptor;
  }
  else{
    a1 = AA;
    a2 = Acceptor;
    a3 =  AA->Cov_neighbor_list[0] == Acceptor ? AA->Cov_neighbor_list[1] : AA->Cov_neighbor_list[0];
  }

  Math3D::Vector3 normal2 = UnitNormal(a1->m_Position,a2->m_Position, a3->m_Position);

	double angle = VectorAngle(normal1,normal2);

  if (angle < Pi/2.0 )
    angle = Pi - angle;

//  cout<<"Out of plane angle is: "<<angle<<endl;
  return angle;
}

double Hbond::computeEnergy() {

  // Mayo energy function, from Dahiyat, Gordon, and Mayo (1997).
  // Automated design of the surface positions of protein helices. Protein Science 6: 1333-1337.

  const double d0 = 8.0; //energy well-depth
  const double r0 = 2.8; //h-bond equilibrium distance
  const double psi0 = toRadian(109.5); // sp3 optimal angle

  double distance_D_A = getDistance_D_A();
  double ratioSquared = r0 * r0 / distance_D_A / distance_D_A;
  double energyDist = d0 * (5 * pow(ratioSquared, 6) - 6 * pow(ratioSquared, 5)); //distance-dependent part
//  log("report")<<"Distance energy part: "<<energyDist<<endl;

  double theta = getAngle_D_H_A();

  //Now we compute the angular term, dependent on the four different cases
  //Todo: CHECK IF GEOMETRIC CRITERIA ARE FULFILLED

  double angularEnergy = cos(theta) * cos(theta); //This is a factor present in all cases
  double energy, psi, phi;

  /// Case 1: donor sp3 and acceptor sp3
  if (m_D_sp3 && m_A_sp3 ) {
//    log("report")<<"Using case D_sp3 A_sp3"<<endl;
    psi = getAngle_H_A_AA();
    energy = energyDist * angularEnergy * cos(psi - psi0) * cos(psi - psi0);
//    log("report")<<"Energy: "<<energy<<", initial energy: "<<m_iniEnergy<<endl;
  }
    /// Case 2: donor sp3 and acceptor sp2
  else if (m_D_sp3 && m_A_sp2 ) {
//    log("report")<<"Using case D_sp3 A_sp2"<<endl;
    psi = getAngle_H_A_AA();
    energy = energyDist * angularEnergy * cos(psi) * cos(psi);
//    log("report")<<"Energy: "<<energy<<", initial energy: "<<m_iniEnergy<<endl;
  }
    /// Case 3: donor sp2 and acceptor sp3
  else if (m_D_sp2 && m_A_sp3 ) {
//    log("report")<<"Using case D_sp2 A_sp3"<<endl;
    energy = energyDist * angularEnergy * angularEnergy;
//    log("report")<<"Energy: "<<energy<<", initial energy: "<<m_iniEnergy<<endl;
  }
    /// Case 4: donor sp2 and acceptor sp2
  else if (m_D_sp2 && m_A_sp2 ) {
//    log("report")<<"Using case D_sp2 A_sp2"<<endl;
    phi = getOutOfPlaneAngle();
    psi = getAngle_H_A_AA();
    psi = max(psi, phi);
    energy = energyDist * angularEnergy * cos(psi) * cos(psi);
//    log("report")<<"Energy: "<<energy<<", initial energy: "<<m_iniEnergy<<endl;
  }

//  log("report")<<"Donor SP2: "<<m_D_sp2<<", Acceptor SP2: "<<m_A_sp2<<", Donor SP3: "<<m_D_sp3<<", Acceptor SP3: "<<m_A_sp3<<endl;

  return energy;
}

void Hbond::identifyHybridization() {

  m_D_sp2 = false;
  m_A_sp2 = false;
  m_D_sp3 = false;
  m_A_sp3 = false;

  /// Using a simple cut-off value of 115 degrees to separate sp2 from sp3 hybridization.
  /// Meng and Lewis (1991). Determination of Molecular Topology and Atomic Hybridization States
  /// from Heavy Atom Coordinates.

  const double cutoff_sp2_sp3 = Math::DtoR(115.0);

  double angleD = 0.0;
  Atom *a1, *a2, *a3;
  //Donor

//  a3 = Donor->Cov_neighbor_list[0] == Hatom ? Donor->Cov_neighbor_list[1] : Donor->Cov_neighbor_list[0];
//  angleD = VectorAngle(Hatom->m_Position - Donor->m_Position, a3->m_Position - Donor->m_Position);

  if(Donor->Cov_neighbor_list.size() >= 3) {
    a2 = Donor->Cov_neighbor_list[0];
    a3 = Donor->Cov_neighbor_list[1];
    if (a2 == Hatom)
      a2 = Donor->Cov_neighbor_list[2];
    if (a3 == Hatom)
      a3 = Donor->Cov_neighbor_list[2];
    a1 = Hatom;

    angleD = Angle(a1->m_Position, Donor->m_Position, a2->m_Position);
    angleD += Angle(a2->m_Position, Donor->m_Position, a3->m_Position);
    angleD += Angle(a3->m_Position, Donor->m_Position, a1->m_Position);

    angleD = angleD/3.0; //mean angle
  }
  else{
    a1 = Donor->Cov_neighbor_list[0] == Hatom ? Donor->Cov_neighbor_list[1] : Donor->Cov_neighbor_list[0];
    angleD = Angle(a1->m_Position, Donor->m_Position, Hatom->m_Position);
  }

  log("report")<<"Donor angle at hbond "<<Hatom->getId()<<", "<<Acceptor->getId()<<": "<<Math::RtoD(angleD)<<endl;
  if(angleD <= cutoff_sp2_sp3 )
    m_D_sp3 = true;
  else
    m_D_sp2 = true;

  //Acceptor
  double angleA = 0.0;

//  a3 = AA->Cov_neighbor_list[0] == Acceptor ? AA->Cov_neighbor_list[1] : AA->Cov_neighbor_list[0];
//  angleA = VectorAngle(Acceptor->m_Position - AA->m_Position, a3->m_Position - AA->m_Position);

  if(AA->Cov_neighbor_list.size() >= 3) {
    a2 = AA->Cov_neighbor_list[0];
    a3 = AA->Cov_neighbor_list[1];
    if (a2 == Acceptor)
      a2 = AA->Cov_neighbor_list[2];
    if (a3 == Acceptor)
      a3 = AA->Cov_neighbor_list[2];
    a1 = Acceptor;

    angleA = Angle(a1->m_Position, AA->m_Position, a2->m_Position);
    angleA += Angle(a2->m_Position, AA->m_Position, a3->m_Position);
    angleA += Angle(a3->m_Position, AA->m_Position, a1->m_Position);

    angleA = angleA/3.0; //mean angle
  }
  else{
    a1 = AA->Cov_neighbor_list[0] == Acceptor ? AA->Cov_neighbor_list[1] : AA->Cov_neighbor_list[0];
    angleA = Angle(a1->m_Position, AA->m_Position, Acceptor->m_Position);
  }

  log("report")<<"Acceptor angle at hbond "<<Hatom->getId()<<", "<<Acceptor->getId()<<": "<<Math::RtoD(angleA)<<endl;
  if(angleA <= cutoff_sp2_sp3 )
    m_A_sp3 = true;
  else
    m_A_sp2 = true;

}

bool Hbond::evaluateGeometry() {

  /// Todo: Implement geometric criteria for boundaries of Mayo energy function
  return true;

}