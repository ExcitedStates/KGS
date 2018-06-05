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


#include "HydrophobicBond.h"
#include "core/Atom.h"
#include "math/MathUtility.h"
#include "Logger.h"

using namespace Math3D;
using namespace std;

HydrophobicBond::HydrophobicBond(Atom* c, Atom* s/*, double energy*/): Bond(c, s, "HYB")
{

  m_atom1=c;
  m_atom2=s;
  m_iniDist = getLength();
  m_bars = 1;

  rigidified=false;

  }
HydrophobicBond::HydrophobicBond(HydrophobicBond & hydrophobicBond) {

c=hydrophobicBond.c;
s=hydrophobicBond.s;


    m_atom1 = hydrophobicBond.m_atom1;
    m_atom2 = hydrophobicBond.m_atom2;
    m_bondType = hydrophobicBond.m_bondType;
    m_bars = hydrophobicBond.m_bars;
    rigidified = hydrophobicBond.rigidified;


}

bool HydrophobicBond::isSame (HydrophobicBond * b2) {
    if ( c->getName() == b2->c->getName() &&
         c->getResidue()->getId() == b2->c->getResidue()->getId() &&
         s->getName()==b2->s->getName() &&
         s->getResidue()->getId() == b2->s->getResidue()->getId()){
        return true;
    }
    return false;
}

double HydrophobicBond::getLength() {
  double length = m_atom1->m_position.distanceTo(m_atom2->m_position);
  return length;
}
