//
// Created by Dominik Budday on 02.08.18.
//

#ifndef KGS_DEEROPTIONS_H
#define KGS_DEEROPTIONS_H

#include <string>
#include <vector>
#include <stdio.h>
#include "Selection.h"
#include "ApplicationOptions.h"

class Molecule;

class Atom;

class DeerOptions : ApplicationOptions {

 public:
  DeerOptions();
  DeerOptions(int argc, char *argv[]);

  /** The working directory */
  std::string workingDirectory;
  /** File-path for the initial structure. */
  std::string initialStructureFile;
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
//  /** Desired norm of step to next sample, can be decreased during collision with decreaseSteps and decreaseFactor */
//  bool flexibleRibose;
  /** Desired metric */
  std::string metric_string;
  /** Selection-pattern passed to metric */
  std::string metricSelection;
  /** Desired planner */
  std::string planner_string;
  /** Percentage to bias random sample to target conf*/
  double biasToTarget;
  /** List of residues supposed active in transition*/
  std::string residueNetwork;
  /** The random generator seed. */
  int seed;
  /** Save output: Indicate the level of files that shall be saved. */
  int saveData;
  /** Percentage to bias random sample to target conf*/
  double convergeDistance;
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
  /** The size of the propagating list of configurations. */
  int frontSize;
  /** Cut-off for svd computation (magnitude of smallest singular value in the nullspace)*/
  double svdCutoff;
  /** Option for collapsing rigid edges. */
  int collapseRigid;
  /** Specified distance to reach between couple of atoms */
  std::string relativeDistances;
  /** Indicates if constraint strain should be predicted and printed */
  bool predictStrain;
  /** How far should samples be explored after DEER distances have been reached. */
  double explorationRadius;
  /** Explore conformations after DEER distances have been reached */
  bool explore;

  void print();

 private:

  void printUsage(char *pname);

  void initializeVariables();
};


#endif //KGS_DEEROPTIONS_H
