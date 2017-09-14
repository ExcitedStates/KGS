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

#include <cassert>
#include <cmath>
#include "KinVertex.h"

#include "Logger.h"

using namespace std;

KinVertex::KinVertex (Rigidbody* rb_ptr):
    m_rigidbody(rb_ptr)
{
  m_parent = nullptr;
  Visited = false;

  if(rb_ptr!=nullptr)
    rb_ptr->setVertex(this);

  m_transformation.setIdentity();
}

KinVertex::~KinVertex () {
  for (auto eit=m_edges.begin(); eit!=m_edges.end(); ++eit) {
    delete *eit;
  }
  if(m_rigidbody)
    m_rigidbody->setVertex(nullptr);
}

void KinVertex::setParent(KinVertex* v) {
  m_parent = v;
}

void KinVertex::addEdge (KinEdge *edge) {
  m_edges.push_back( edge );
}

KinEdge* KinVertex::findEdge(const KinVertex* v) const
{
  for(auto const& edge: m_edges){
    if( edge->EndVertex==v )
      return edge;
  }
  return nullptr;
}


void KinVertex::print() const {
  log() << "KinVertex";
  if(m_rigidbody)
    for (vector<Atom*>::iterator it=m_rigidbody->Atoms.begin(); it!=m_rigidbody->Atoms.end(); ++it)
      log() << (*it)->getId() << "+";
  log() << endl;
}

void KinVertex::forwardPropagate()
{
  for(auto const& edge: m_edges){
    edge->forwardPropagate();
  }

  //Apply transformation AFTER propagation. This is important.
  transformAtoms();
}

void KinVertex::transformAtoms()
{
  if(m_rigidbody==nullptr) return;

  for (auto const& atom: m_rigidbody->Atoms){
    Math3D::Vector3 newPos = m_transformation * atom->m_referencePosition;

    atom->m_position.x = newPos.x;
    atom->m_position.y = newPos.y;
    atom->m_position.z = newPos.z;

    assert( !std::isnan(newPos.x) );
    assert( !std::isnan(newPos.y) );
    assert( !std::isnan(newPos.z) );
  }
}

