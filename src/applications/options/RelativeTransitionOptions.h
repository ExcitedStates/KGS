#ifndef RELATIVETRANSITIONOPTIONS_H
#define RELATIVETRANSITIONOPTIONS_H

#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"

class RelativeTransitionOptions
{
 public:

  /** File-path for the initial structure. */
  std::string initialStructureFile;
  /** The working directory */
  std::string workingDirectory;
  /** File containing hydrogen bond definitions. */
  std::string hydrogenbondFile;
  /** Method of identifying hydrogen bonds. */
  std::string hydrogenbondMethod;
  /** Atom-pairs that should be treated as h-bonds. */
  std::vector<std::string> extraCovBonds;

  /** How many samples should be generated */
  int samplesToGenerate;
  /** Defines the method used to calculate a perturbation or gradient */
  int gradient;
  /** Collision factor indicating the multiplication factor of all atoms van der Waals radii. */
  double collisionFactor;
  /** If a non-colliding structure can not be found, try this many times to decrease the stepsize */
  int decreaseSteps;
  /** If a non-colliding structure can not be found, decrease the stepsize by this factor */
  double decreaseFactor;
  /** The largest allowable change in torsion angle */
  double maxRotation;
  /** step size factor*/
  double stepSize;
  /** Desired metric */
  std::string metric_string;
  /** Selection-pattern passed to metric */
  std::string metricSelection;

  /** The random generator seed. */
  int seed;
  /** Percentage to bias random sample to target conf*/
  double biasToTarget;
  /** List of residues supposed active in transition*/
  std::string residueNetwork;
  /** Percentage to bias random sample to target conf*/
  double convergeDistance;
  /** Uses clashing atoms to define additional constraints. */
  std::string selectAtoms;
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
  /** Cut-off for svd computation (magnitude of smallest singular value in the nullspace)*/
  double svdCutoff;
  /** Option for collapsing rigid edges. */
  int collapseRigid;
  /** Specified distance to reach between couple of atoms */
  std::string relativeDistances;

  void print();

 private:
  RelativeTransitionOptions();
  RelativeTransitionOptions(int argc, char* argv[] );

  void printUsage(char* pname);
  inline bool fileExists ( const std::string& name);
  void initializeVariables();

  //Singleton pattern
  static RelativeTransitionOptions* instance;
 public:
  static RelativeTransitionOptions* getOptions();
  static RelativeTransitionOptions* createOptions(int argc, char* argv[] );
  static RelativeTransitionOptions* createOptions();

};

#endif // SAMPLINGOPTIONS_H
