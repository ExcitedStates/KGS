/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

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


#include "HBond.h"
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
	m_atom1 = hatom;
	m_atom2 = acceptor;
	m_bondType = "HB";
	m_bars = 5;
	rigidified = false;
	m_iniDist_H_A = VectorLength(m_atom1->m_position,m_atom2->m_position);
	m_iniAngle_D_H_A = getAngle_D_H_A();
	m_iniAngle_H_A_AA = getAngle_H_A_AA();

  identifyHybridization();

  m_iniEnergy = DEFAULT_HBOND_ENERGY;
  m_iniEnergy = computeEnergy();
}

Hbond::Hbond(Hbond & hbond) {
	Hatom = hbond.Hatom;
	Acceptor = hbond.Acceptor;
	Donor = hbond.Donor;
	AA = hbond.AA;
	m_iniEnergy = hbond.m_iniEnergy;

	// super-class attributes
	m_atom1 = hbond.m_atom1;
	m_atom2 = hbond.m_atom2;
	m_bondType = hbond.m_bondType;
	m_bars = hbond.m_bars;
	rigidified = hbond.rigidified;
	m_iniDist_H_A = hbond.m_iniDist_H_A;
	m_iniAngle_D_H_A = hbond.m_iniAngle_D_H_A;
	m_iniAngle_H_A_AA = hbond.m_iniAngle_H_A_AA;

  m_acceptorHybridization = hbond.m_acceptorHybridization;
  m_donorHybridization = hbond.m_donorHybridization;

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
//    return Acceptor->m_position + (x*idealH[0]) + (y*idealH[1]) + (z*idealH[2]);
//}
//
//Vector3 Hbond::getIdealAcceptorPoint(){
//    Vector3 x,y,z;
//    coordinateSystem(Hatom, x,y,z);
//    return Hatom->m_position + (x*idealA[0]) + (y*idealA[1]) + (z*idealA[2]);
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
//    x = a1->m_position -a2->m_position;
//    y = a3->m_position -a2->m_position;
//    z = cross(x,y);
//}

double Hbond::getLength() {
  double length = Hatom->m_position.distanceTo(Acceptor->m_position);
  return length;
}

double Hbond::getDistance_D_A() {
  double distance = Donor->m_position.distanceTo(Acceptor->m_position);
  return distance;
}

double Hbond::getAngle_D_H_A() {
  double ret = Angle(Donor->m_position, Hatom->m_position, Acceptor->m_position);
  return ret;
}

double Hbond::getAngle_H_A_AA() {
  double ret = Angle(Hatom->m_position, Acceptor->m_position, AA->m_position);
  return ret;
}

double Hbond::getOutOfPlaneAngle() {

	//Todo: Check if one should rather use the other covalent neighbors instead of H-atom and Acceptor

	Atom *a1, *a2, *a3;

//  if(Donor->Cov_neighbor_list.size() >= 3) {
//    a1 = Donor->Cov_neighbor_list.at(0);
//    a2 = Donor->Cov_neighbor_list.at(1);
//    a3 = Donor->Cov_neighbor_list.at(2);
//  }
//  else
  if(Donor->Cov_neighbor_list.size() >= 2){
    a1 = Donor;
    a2 = Donor->Cov_neighbor_list.at(0);
    a3 = Donor->Cov_neighbor_list.at(1);
  }
  else{
    a1 = Donor;
    a2 = Hatom;
    a3 = Acceptor;
  }

  Math3D::Vector3 normal1 = UnitNormal(a1->m_position,a2->m_position, a3->m_position);

//  if(Acceptor->Cov_neighbor_list.size() >= 3) {
//    a1 = Acceptor->Cov_neighbor_list.at(0);
//    a2 = Acceptor->Cov_neighbor_list.at(1);
//    a3 = Acceptor->Cov_neighbor_list.at(2);
//  }
//  else
  if(AA->Cov_neighbor_list.size() >= 2){
    a1 = AA;
    a2 = AA->Cov_neighbor_list.at(0);
    a3 = AA->Cov_neighbor_list.at(1);
  }
  else{
    a1 = AA;
    a2 = Acceptor;
    a3 = Hatom;
  }

  Math3D::Vector3 normal2 = UnitNormal(a1->m_position,a2->m_position, a3->m_position);

	double angle = VectorAngle(normal1,normal2);

  if (angle < Pi/2.0 )//Return the angle > 90 Degrees
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
  if (m_donorHybridization == 3 && m_acceptorHybridization == 3 ) {
    log("report") << "Hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": Using case D_sp3 A_sp3"<<endl;
    psi = getAngle_H_A_AA();
    energy = energyDist * angularEnergy * cos(psi - psi0) * cos(psi - psi0);
  }
    /// Case 2: donor sp3 and acceptor sp2
  else if (m_donorHybridization == 3 && m_acceptorHybridization == 2 ) {
    log("report") << "Hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": Using case D_sp3 A_sp2"<<endl;
    psi = getAngle_H_A_AA();
    energy = energyDist * angularEnergy * cos(psi) * cos(psi);
  }
    /// Case 3: donor sp2 and acceptor sp3
  else if (m_donorHybridization == 2 && m_acceptorHybridization == 3 ) {
    log("report") << "Hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": Using case D_sp2 A_sp3"<<endl;
    energy = energyDist * angularEnergy * angularEnergy;
  }
    /// Case 4: donor sp2 and acceptor sp2
  else if (m_donorHybridization == 2 && m_acceptorHybridization == 2 ) {
    log("report") << "Hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": Using case D_sp2 A_sp2"<<endl;
    phi = getOutOfPlaneAngle();
    psi = getAngle_H_A_AA();
    psi = max(psi, phi);
    energy = energyDist * angularEnergy * cos(psi) * cos(psi);
  }
  if(m_iniEnergy == DEFAULT_HBOND_ENERGY){
    log("report")<<"Energy: "<<energy<<", initial energy: "<<energy<<endl;
  }
  else{
    log("report")<<"Energy: "<<energy<<", initial energy: "<<m_iniEnergy<<endl;
  }

  log("report")<<"Donor SP: "<<m_donorHybridization<<", Acceptor SP: "<<m_acceptorHybridization<<endl;
  return energy;
}

void Hbond::identifyHybridization() {

  /// ToDo: Change bool values to integers, either 2 or 3
  /// For the acceptor hybridization include the hydrogen as well, because
  /// it also binds an electron
  /// Then we should be good to go!

  /// Using a simple cut-off value of 115 degrees to separate sp2 from sp3 hybridization.
  /// Meng and Lewis (1991). Determination of Molecular Topology and Atomic Hybridization States
  /// from Heavy Atom Coordinates.

  const double cutoff_sp2_sp3 = Math::DtoR(115.0);
  const double reviewAngle = Math::DtoR(122.0);
  const double dist_C_sp2_sp3 = 1.43;
  const double dist_N_sp2_sp3 = 1.40;

  double angleD = 0.0;
  Atom *a1, *a2, *a3;
  //Donor
  log("report")<<endl;
//  int donorNeighbors = Donor->Cov_neighbor_list.size();
  std::vector<Atom*> d_consideredNeighbors = Donor->Cov_neighbor_list;//Donor->heavyAtomNeighbors();
  int donorNeighbors = d_consideredNeighbors.size();

  log("report") << "Hybridisation Donor for hbond: " << Hatom->getId() << ", " << Acceptor->getId() << endl;
  log("report") << "HAV: " << donorNeighbors <<", # neighbors: "<<Donor->Cov_neighbor_list.size()<<endl;

  if (donorNeighbors > 3) {
    m_donorHybridization = 3;//Don't need to check for angles anymore, it's sp3
    log("report") << "Donor is sp3, no evaluation. "<< endl;
  }
  else if(donorNeighbors < 2){//Does this hold for oxygen as well
    m_donorHybridization = 2;//Don't need to check for angles anymore, it's sp2
    log("report") << "Donor is sp2, no evaluation"<< endl;
  }
  else {
    if (donorNeighbors == 3) {
      a1 = d_consideredNeighbors[0];
      a2 = d_consideredNeighbors[1];
      a3 = d_consideredNeighbors[2];

      angleD = Angle(a1->m_position, Donor->m_position, a2->m_position);
      angleD += Angle(a2->m_position, Donor->m_position, a3->m_position);
      angleD += Angle(a3->m_position, Donor->m_position, a1->m_position);

      angleD = angleD / 3.0; //mean angle

      if (angleD <= cutoff_sp2_sp3) //Cut-off angle from Chimera
        m_donorHybridization = 3;
      else
        m_donorHybridization = 2;
      log("report") << "Avg angle: "<<Math::RtoD(angleD) << endl;
      log("report") << "Donor is sp "<<m_donorHybridization<<endl;
    }
    else {//(donorNeighbors == 2)
      if(Donor->getType() == "O" || Donor->getType() == "S"){//these must be sp3
        m_donorHybridization = 3;
        log("report") << "Donor is sp3, O or S atom. "<< endl;
      }
      else{//treatment of N (and C)
        a1 = d_consideredNeighbors[0];
        a2 = d_consideredNeighbors[1];
        angleD = Angle(a1->m_position, Donor->m_position, a2->m_position);
        log("report") << "Single donor angle at hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": "
                      << Math::RtoD(angleD) << endl;
        if (angleD > reviewAngle) { //Cut-off angle from Chimera for further consideration
          m_donorHybridization = 2;
          log("report") << "Donor is sp2, review angle passed. "<< endl;
        }
        else{//evaluate further criteria to check if sp2 or sp3
          log("report") << "Donor under distance review. "<< endl;
          double distance1 = Donor->distanceTo(d_consideredNeighbors[0]);
          double distance2 = Donor->distanceTo(d_consideredNeighbors[1]);
          double avgDistance = (distance1+distance2)/2;
          log("report") << "Avg distance is "<<avgDistance<< endl;

          if(Donor->getType() == "C"){
            if(avgDistance < dist_C_sp2_sp3)
              m_donorHybridization = 2;
            else
              m_donorHybridization = 3;
          }
          if(Donor->getType() == "N"){
            if(avgDistance < dist_N_sp2_sp3)
              m_donorHybridization = 2;
            else
              m_donorHybridization = 3;
          }
          log("report") << "Donor is sp "<<m_donorHybridization<<endl;
        }
      }
    }
  }
//  log("report") << "Donor angle at hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": "
//                << Math::RtoD(angleD) << endl;
//  log("report") << "Donor is sp "<<m_donorHybridization<<endl;

  //Acceptor
  double angleA = 0.0;

//  int acceptorNeighbors = Acceptor->Cov_neighbor_list.size();
  std::vector<Atom*> a_consideredNeighbors = Acceptor->Cov_neighbor_list;//Acceptor->heavyAtomNeighbors();
  int acceptorNeighbors = a_consideredNeighbors.size();

  log("report") << "Hybridisation Acceptor for hbond: " << Hatom->getId() << ", " << Acceptor->getId() << endl;
  log("report") << "HAV: " << acceptorNeighbors <<", # neighbors: "<<Acceptor->Cov_neighbor_list.size()<<endl;

  //The acceptor is non-covalently bound to the hydrogen, which needs to be included
  //in calculating the hybridization state, because they share electrons via the h-bond.
  //For the acceptor, we include the non-covalently bound hydrogen atom in the count,
  //as this influences desired angles of the h-bond
//  acceptorNeighbors += 1;

  if (acceptorNeighbors > 3) {
    m_acceptorHybridization = 3;//Don't need to check for angles anymore, it's sp3
    log("report") << "Acceptor is sp3, no evaluation. "<< endl;
  }
  else if(acceptorNeighbors < 2){//Does this hold for oxygen as well
    m_acceptorHybridization = 2;//Don't need to check for angles anymore, it's sp2
    log("report") << "Acceptor is sp2, no evaluation"<< endl;
  }
  else {
    if (acceptorNeighbors == 3) {
      a1 = Acceptor->Cov_neighbor_list[0];
      a2 = Acceptor->Cov_neighbor_list[1];
      a3 = Acceptor->Cov_neighbor_list[2];

      angleA = Angle(a1->m_position, Acceptor->m_position, a2->m_position);
      angleA += Angle(a2->m_position, Acceptor->m_position, a3->m_position);
      angleA += Angle(a3->m_position, Acceptor->m_position, a1->m_position);

      angleA = angleA / 3.0; //mean angle

      if (angleA <= cutoff_sp2_sp3) //Cut-off angle from Chimera
        m_acceptorHybridization = 3;
      else
        m_acceptorHybridization = 2;
      log("report") << "Avg angle: " << Math::RtoD(angleD) << endl;
      log("report") << "Acceptor is sp " << m_donorHybridization << endl;
    }
    else {//(acceptorNeighbors == 2)
      if(Acceptor->getType() == "O" || Acceptor->getType() == "S"){//these must be sp3
        m_acceptorHybridization = 3;
        log("report") << "Donor is sp3, O or S atom. "<< endl;
      }
      else {//treatment of N (and C)
        a1 = a_consideredNeighbors[0];
        a2 = a_consideredNeighbors[1];
        angleA = Angle(a1->m_position, Acceptor->m_position, a2->m_position);
        log("report") << "Single acceptor angle at hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": "
                      << Math::RtoD(angleA) << endl;
        if (angleA > reviewAngle) { //Cut-off angle from Chimera for further consideration
          m_acceptorHybridization = 2;
          log("report") << "Acceptor is sp2, review angle passed. " << endl;
        } else {//evaluate further criteria to check if sp2 or sp3
          log("report") << "Acceptor under distance review. " << endl;
          double distance1 = Acceptor->distanceTo(a_consideredNeighbors[0]);
          double distance2 = Acceptor->distanceTo(a_consideredNeighbors[1]);
          double avgDistance = (distance1 + distance2) / 2;
          log("report") << "Avg distance is " << avgDistance << endl;

          if (Acceptor->getType() == "C") {
            if (avgDistance < dist_C_sp2_sp3)
              m_acceptorHybridization = 2;
            else
              m_acceptorHybridization = 3;
          }
          if (Acceptor->getType() == "N") {
            if (avgDistance < dist_N_sp2_sp3)
              m_acceptorHybridization = 2;
            else
              m_acceptorHybridization = 3;
          }
          log("report") << "Acceptor is sp " << m_donorHybridization << endl;
        }
      }
    }
  }
    //  log("report") << "Acceptor angle at hbond " << Hatom->getId() << ", " << Acceptor->getId() << ": "
//                << Math::RtoD(angleD) << endl;
//  log("report") << "Acceptor is sp "<<m_acceptorHybridization<<endl;
  log("report")<<endl;
}

bool Hbond::evaluateGeometry() {

  /// Todo: Implement geometric criteria for boundaries of Mayo energy function
  return true;

}
