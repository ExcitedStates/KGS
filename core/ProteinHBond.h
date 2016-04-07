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
  double Energy;
  double Dist_H_A;
  double Ang_H_A_AA;
  double iniLength;
  double iniOrientLeft;
  double iniOrientRight;

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
  double getLeftAngle();
  double getRightAngle();
 private:
  void coordinateSystem(Atom* a, Math3D::Vector3& x, Math3D::Vector3& y, Math3D::Vector3& z );
  Math3D::Vector3 idealH;
  Math3D::Vector3 idealA;

};

#endif
