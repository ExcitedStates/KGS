//
// Created by Dominik Budday on 19.12.16.
//

#ifndef KGS_RIGIDITYOPTIONS_H
#define KGS_RIGIDITYOPTIONS_H



#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"

class Molecule;

class Atom;

class RigidityOptions {
 public:

  /** The working directory */
  std::string workingDirectory;
  /** File-path for the initial structure. */
  std::string initialStructureFile;
  /** Annotation file-path common for all initial structures. */
  std::string annotationFile;
  /** File containing hydrogen bond definitions. */
  std::string hydrogenbondFile;
  /** Method of identifying hydrogen bonds. */
  std::string hydrogenbondMethod;
  /** Atom-pairs that should be treated as h-bonds. */
  std::vector<std::string> extraCovBonds;
  /** Collision factor indicating the multiplication factor of all atoms van
        der Waals radii. */
  double collisionFactor;
  /** The random generator seed. */
  int seed;
  /** Save output: Indicate the level of files that shall be saved. */
  int saveData;
  /** List of residues supposed active in transition, TODO: allow moveable subset in rigidity analysis as well*/
  std::string residueNetwork;
  /** Atom selection for gradients, distance computation etc., e.g. heavy, name Ca, backbone, all. Default is heavy*/
  bool preventClashes;
  /** The vector of atom IDs that will be part of root rigid bodies */
  std::vector<int> roots;
  /** Atoms used for collision detection and clash constraints: "all, backbone, heavy" */
  std::string collisionCheck;
  /** Cut-off for svd computation (magnitude of smallest singular value in the nullspace)*/
  double svdCutoff;
  /** Option for collapsing rigid edges. */
  int collapseRigid;

  void print();

 private:
  RigidityOptions();

  RigidityOptions(int argc, char *argv[]);

  void printUsage(char *pname);

  inline bool fileExists(const std::string &name);

  void initializeVariables();

  //Singleton pattern
  static RigidityOptions *instance;
 public:
  static RigidityOptions *getOptions();

  static RigidityOptions *createOptions(int argc, char *argv[]);

  static RigidityOptions *createOptions();

};
#endif //KGS_RIGIDITYOPTIONS_H
