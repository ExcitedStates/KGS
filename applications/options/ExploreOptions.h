#ifndef EXPLOREOPTIONS_H
#define EXPLOREOPTIONS_H

#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"

class Molecule;

class Atom;

class ExploreOptions {
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

  /** How many samples should be generated */
  int samplesToGenerate;
  /** How far from the center structure should samples be explored. */
  double explorationRadius;
  /** Generate new samples from randomly chosen seed samples (instead from last accepted sample). */
  bool sampleRandom;
  /** Defines the method used to calculate a perturbation or gradient */
  int gradient;
  /** Collision factor indicating the multiplication factor of all atoms van
        der Waals radii. */
  double collisionFactor;
  /** If a non-colliding structure can not be found, try this many times to decrease the stepsize */
  int decreaseSteps;
  /** If a non-colliding structure can not be found, decrease the stepsize by this factor */
  double decreaseFactor;
  /** The largest allowable change in torsion angle */
  double maxRotation;
  /** step size factor*/
  double stepSize;
//  /** Desired norm of step to next sample, can be decreased during collision with decreaseSteps and decreaseFactor */
//  bool flexibleRibose;
//  /** Residue selection of residues to rebuild with exactIK */
//  std::string exactIKselection;
  /** Desired metric */
  std::string metric_string;
  /** Selection-pattern passed to metric */
  std::string metricSelection;
  /** Desired planner */
  std::string planner_string;
  /** Generate new samples from randomly chosen seed samples (instead from last accepted sample). */
  bool scaleToRadius;
//  /** Percentage to bias random sample to target conf*/
//  double biasToTarget;
  /** Number perturbations that the Poisson planner tries before a sample is 'closed' */
  int poissonMaxRejectsBeforeClose;

//  /** How frequently to perform a rebuild move. Should be a number betwen 0 and 1. */
//  double rebuild_frequency;
//
//  /** The number of initial structures to add by iteratively choosing a random initial structure and
//    rebuilding all its free loop-regions. */
//  int rebuildInitialStructures;
//
//  /** The aggression of the rebuild procedure. */
//  int rebuildAggression;

  /** The random generator seed. */
  int seed;
  /** Save output: Indicate the level of files that shall be saved. */
  int saveData;
//  /** Sample in reverse direction as well */
//  bool sampleReverse;
  /** List of residues supposed active in transition*/
//  std::vector<int> residueNetwork;
  std::string residueNetwork;
//  /** Percentage to bias random sample to target conf*/
//  double convergeDistance;
  /** Align configs to initial. */
  bool alignAlways;
  /** Align initial and target configuration in the beginning. */
  bool alignIni;
  /** Atom selection for gradients, distance computation etc., e.g. heavy, name Ca, backbone, all. Default is heavy*/
  bool preventClashes;
  /** Specifies the residues of the molecule that will undergo RMSD alignment during sampling. */
  std::string alignSelection;
  /** Specifies the residues used for gradient computation. */
  std::string gradientSelection;
  /** The vector of atom IDs that will be part of root rigid bodies */
  std::vector<int> roots;
  /** Whether or not to project the gradient onto the constraint manifold. */
  bool projectConstraints;
  /** Atoms used for collision detection and clash constraints: "all, backbone, heavy" */
  std::string collisionCheck;
//  /** The size of the propagating list of configurations. */
//  int frontSize;
//  /** Max number of trials before switching search directions. */
//  int switchAfter;
  /** Cut-off for svd computation (magnitude of smallest singular value in the nullspace)*/
  double svdCutoff;
  /** Option for collapsing rigid edges. */
  int collapseRigid;

//  /** Specified distance to reach between couple of atoms */
//  std::string relativeDistances;

  void print();

//  void setResidueNetwork(const Molecule * protein);
//  void setAtomSets(const Molecule * protein, Molecule * target);
//  const std::vector<Atom*>* getAtomsAlign() const { return &m_atomsAlign;} //pointer return, as nullptr required
//  const std::vector<Atom*>* getAtomsMoving() const { return &m_atomsMoving;} //pointer return, as nullptr required

 private:
  ExploreOptions();

  ExploreOptions(int argc, char *argv[]);

  void printUsage(char *pname);

  inline bool fileExists(const std::string &name);

  void initializeVariables();

//  std::vector<Atom*> m_atomsAlign; //atoms used for alignment, specified via alignSelection
//  std::vector<Atom*> m_atomsMoving; //atoms used for gradient computation, specified via gradientSelection

  //Singleton pattern
  static ExploreOptions *instance;
 public:
  static ExploreOptions *getOptions();

  static ExploreOptions *createOptions(int argc, char *argv[]);

  static ExploreOptions *createOptions();

};

#endif // SAMPLINGOPTIONS_H
