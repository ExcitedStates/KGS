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
#ifndef HBOND_H
#define HBOND_H

#include "core/Bond.h"
#include "math3d/primitives.h"

#define DEFAULT_HBOND_ENERGY -100
#define DEFAULT_HBOND_DIST 100
#define DEFAULT_HBOND_ANG 100

class Hbond : public Bond {
 public:

  Atom* Donor;
  Atom* Hatom;
  Atom* Acceptor;
  Atom* AA;

  Hbond(Atom* hatom, Atom* acceptor, Atom* donor, Atom* aa, double energy=DEFAULT_HBOND_ENERGY);
  Hbond(Hbond & hbond);
  bool isSame (Hbond * b2);
  Atom* atom1(){ return Bond::Atom1; }
  Atom* atom2(){ return Bond::Atom2; }
  Math3D::Vector3 getIdealHPoint();
  Math3D::Vector3 getIdealAcceptorPoint();
  double getLength();
  double getDistance_D_A();
  double getAngle_D_H_A();
  double getAngle_H_A_AA();
  double getOutOfPlaneAngle();
  bool evaluateGeometry();
  double computeEnergy();

  void identifyHybridization();

// private:
//  void coordinateSystem(Atom* a, Math3D::Vector3& x, Math3D::Vector3& y, Math3D::Vector3& z );
//  Math3D::Vector3 idealH;
//  Math3D::Vector3 idealA;

 private:
//These are initial distance and orientations, only set in the beginning
  double m_iniDist_H_A;
  double m_iniAngle_D_H_A;
  double m_iniAngle_H_A_AA;
  double m_iniEnergy;
  bool m_D_sp2;
  bool m_D_sp3;
  bool m_A_sp2;
  bool m_A_sp3;

  //Getters and Setters for initially determined values
 public:
  double getIniLength() const {
    return m_iniDist_H_A;
  }
  void setIniLength(double iniDist_H_A) {
    m_iniDist_H_A = iniDist_H_A;
  }

  double getIniAngle_D_H_A() const {
    return m_iniAngle_D_H_A;
  }

  void setIniAngle_D_H_A(double iniAngle_D_H_A) {
    m_iniAngle_D_H_A = iniAngle_D_H_A;
  }

  double getIniAngle_H_A_AA() const {
    return m_iniAngle_H_A_AA;
  }

  void setIniAngle_H_A_AA(double iniAngle_H_A_AA) {
    m_iniAngle_H_A_AA = iniAngle_H_A_AA;
  }

  double getIniEnergy() const {
    return m_iniEnergy;
  }

  void setIniEnergy(double energy) {
    m_iniEnergy = energy;
  }

  bool donorSP2() const {
    return m_D_sp2;
  }

  bool donorSP3() const {
    return m_D_sp3;
  }

  bool acceptorSP2() const {
    return m_A_sp2;
  }

  bool acceptorSP3() const {
    return m_A_sp3;
  }

};

#endif
