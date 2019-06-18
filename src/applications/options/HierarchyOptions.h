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
//
// Created by Dominik Budday on 20.12.16.
//

#ifndef KGS_HIERARCHYOPTIONS_H
#define KGS_HIERARCHYOPTIONS_H



#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"

class Molecule;

class Atom;

class HierarchyOptions {
 public:

  /** The working directory */
  std::string workingDirectory;
  /** File-path for the initial structure. */
  std::string initialStructureFile;
  /** Annotation file-path common for all initial structures. */
  std::string targetStructureFile;
  /** Annotation file-path common for all initial structures. */
  std::string annotationFile;
  /** File containing hydrogen bond definitions. */
  std::string hydrogenbondFile;
  /** Method of identifying hydrogen bonds. */
  std::string hydrogenbondMethod;
  /** Atom-pairs that should be treated as h-bonds. */
  std::vector<std::string> extraCovBonds;

  /** Which samples should be generated. Most of the analysis can be done without samples (0 (default) up to dof size)*/
  int samples;
  /** Collision factor indicating the multiplication factor of all atoms van
        der Waals radii. */
  double collisionFactor;
  /** The random generator seed. */
  int seed;
  /** step size factor*/
  double stepSize;
  /** Save output: Indicate the level of files that shall be saved. */
  int saveData;
  /** List of residues active in analysis*/
  std::string residueNetwork;
  /** The vector of atom IDs that will be part of root rigid bodies */
  std::vector<int> roots;
  /** Atoms used for collision detection and clash constraints: "all, backbone, heavy" */
  std::string collisionCheck;
  /** Cut-off for svd computation (magnitude of smallest singular value in the nullspace)*/
  double svdCutoff;
  /** Sink/Source to identify transferred DoF between different areas. */
  std::string sink;
  /** Sink/Source to identify transferred DoF between different areas. */
  std::string source;
  /** Sample lowest free-energy motions (instead of along IDs in the V matrix). */
  bool sampleFree;
  void print();

 private:
  HierarchyOptions();

  HierarchyOptions(int argc, char *argv[]);

  void printUsage(char *pname);

  inline bool fileExists(const std::string &name);

  void initializeVariables();

  //Singleton pattern
  static HierarchyOptions *instance;
 public:
  static HierarchyOptions *getOptions();

  static HierarchyOptions *createOptions(int argc, char *argv[]);

  static HierarchyOptions *createOptions();

};
#endif //KGS_HIERARCHYOPTIONS_H
