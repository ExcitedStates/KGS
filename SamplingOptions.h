#ifndef SAMPLINGOPTIONS_H
#define SAMPLINGOPTIONS_H

#include <string>
#include <vector>
#include <stdio.h>

class Molecule;
class Atom;

class SamplingOptions
{
 public:

  /** The working directory */
  std::string workingDirectory;
  /** File-path for the initial structure. */
  std::string initialStructureFile;
  /** File-path for all initial structures. */
  std::vector<std::string> initialStructureFiles;
  std::string targetStructureFile;
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
  /** Desired norm of step to next sample, can be decreased during collision with decreaseSteps and decreaseFactor */
  bool flexibleRibose;
  /** Desired m_metric */
  std::string metric_string;
  /** Desired planner */
  std::string planner_string;
  /** Generate new samples from randomly chosen seed samples (instead from last accepted sample). */
  bool scaleToRadius;
  /** Percentage to bias random sample to target conf*/
  double biasToTarget;
  /** Number perturbations that the Poisson planner tries before a sample is 'closed' */
  int poisson_max_rejects_before_close;

  /** The length of free (i.e. not annotated as fixed) fragments being rebuilt */
  int rebuild_fragment_length;
  /** How frequently to perform a rebuild move. Should be a number betwen 0 and 1. */
  double rebuild_frequency;

  /** The number of initial structures to add by iteratively choosing a random initial structure and
    rebuilding all its free loop-regions. */
  int rebuildInitialStructures;

  /** The aggression of the rebuild procedure. */
  int rebuildAggression;

  /** The random generator seed. */
  int seed;
  /** Save output: Indicate the level of files that shall be saved. */
  int saveData;
  /** Sample in reverse direction as well */
  bool sampleReverse;
  /** List of residues supposed active in transition*/
  std::vector<int> residueNetwork;
  /** Percentage to bias random sample to target conf*/
  double convergeDistance;
  /** Align configs to initial. */
  bool alignAlways;
  /** Align initial and target configuration in the beginning. */
  bool alignIni;
  /** Uses clashing atoms to define additional constraints. */
  std::string selectAtoms;
  /** Atom selection for gradients, distance computation etc., e.g. heavy, name Ca, backbone, all. Default is heavy*/
  bool preventClashes;
  /** Specifies the residues of the molecule that will undergo RMSD alignment during sampling. */
  std::string selectionAlign;
  /** Specifies the residues of the molecule that will undergo RMSD alignment during sampling. */
  std::string selectionMoving;
  /** The m_root rigid body id. */
  int root;
  /** Whether or not to project the gradient onto the constraint manifold. */
  bool projectConstraints;
  /** Atoms used for collision detection and clash constraints: "all, backbone, heavy" */
  std::string collisionCheck;
  /** The size of the propagating list of configurations. */
  int frontSize;
  /** Max number of trials before switching search directions. */
  int switchAfter;

  void print();

  void setResidueNetwork(const Molecule * protein);
  void setAtomSets(const Molecule * protein, Molecule * target);
  const std::vector<Atom*>* getAtomsAlign() const { return &m_atomsAlign;} //pointer return, as nullptr required
  const std::vector<Atom*>* getAtomsMoving() const { return &m_atomsMoving;} //pointer return, as nullptr required

 private:
  SamplingOptions();
  SamplingOptions(int argc, char* argv[] );

  void printUsage(char* pname);
  inline bool fileExists ( const std::string& name);
  void initializeVariables();

  std::vector<Atom*> m_atomsAlign; //atoms used for alignment, specified via selectionAlign
  std::vector<Atom*> m_atomsMoving; //atoms used for gradient and RMSD computation, specified via selectionMoving

  //Singleton pattern
  static SamplingOptions* instance;
 public:
  static SamplingOptions* getOptions();
  static SamplingOptions* createOptions(int argc, char* argv[] );
  static SamplingOptions* createOptions();

};

#endif // SAMPLINGOPTIONS_H
