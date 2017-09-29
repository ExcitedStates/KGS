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
// Disjoint Set Data Structure
// Author: Emil Stefanov
// Date: 03/28/06
// Implementaton is as described in http://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include <vector>
#include <stddef.h>

class DisjointSets
{
 public:

  // Create an empty DisjointSets data structure
  DisjointSets();
  // Create a DisjointSets data structure with a specified number of elements (with element id's from 0 to count-1)
  DisjointSets(int count);
  // Copy constructor
  DisjointSets(const DisjointSets & s);
  // Destructor
  ~DisjointSets();

  // Find the set identifier that an element currently belongs to.
  // Note: some internal data is modified for optimization even though this method is constant.
  int FindSet(int element) const;
  // Combine two sets (represented by elements) into one. All elements in those two sets will share the same set id that can be gotten using FindSet.
  void Union(int elemId1, int elemId2);
  // Add a specified number of elements to the DisjointSets data structure. The element id's of the new elements are numbered
  // consecutively starting with the first never-before-used elementId.
  void AddElements(int numToAdd);
  // Returns the number of elements currently in the DisjointSets data structure.
  int NumElements() const;
  // Returns the number of sets currently in the DisjointSets data structure.
  int NumSets() const;

 private:

  // Internal Node data structure used for representing an element
  struct Node
  {
    int rank; // This roughly represent the max height of the node in its subtree
    int index; // The index of the element the node represents
    Node* parent; // The m_parent node of the node
  };

  int m_numElements; // the number of elements currently in the DisjointSets data structure.
  int m_numSets; // the number of sets currently in the DisjointSets data structure.
  std::vector<Node*> m_nodes; // the list of nodes representing the elements
};
