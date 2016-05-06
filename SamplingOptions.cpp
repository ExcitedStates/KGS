#include "SamplingOptions.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "Util.h"
#include "Logger.h"
#include "Selection.h"
#include "core/Molecule.h"
#include "core/Atom.h"
#include "core/Chain.h"

using namespace std;

SamplingOptions::SamplingOptions(){
	initializeVariables();
}

SamplingOptions::SamplingOptions(int argc, char* argv[]):
        extraCovBonds()
{
	initializeVariables();

	if(argc<2){
		enableLogger("so");
		printUsage(argv[0]);
		exit(-1);
	}

	int i;
	for(i=1;i<argc;i++){
		string arg = argv[i];
		if(arg=="--initial"){                       initialStructureFile = argv[++i];                   continue; }
		if(arg=="--target"){                        targetStructureFile = argv[++i];                    continue; }
		if(arg=="--annotation"){                    annotationFile = argv[++i];                         continue; }
		if(arg=="--hbondMethod"){                   hydrogenbondMethod = argv[++i];                     continue; }
		if(arg=="--hbondFile"){                     hydrogenbondFile = argv[++i];                       continue; }
		if(arg=="--extraCovBonds"){                 Selection::split( argv[++i],",", extraCovBonds );   continue; }
		if(arg=="--workingDirectory"){              workingDirectory = argv[++i];                       continue; }
		//if(arg=="--rigidbodiesFromFIRST"){          hydrogenbondFile = argv[++i];                       continue; }
		if(arg=="--samples" || arg=="-s"){          samplesToGenerate = atoi(argv[++i]);                continue; }
		if(arg=="--radius" || arg=="-r"){           explorationRadius = atof(argv[++i]);                continue; }
		if(arg=="--sampleRandom"){                  sampleRandom = Util::stob(argv[++i]);               continue; }
		if(arg=="--scaleToRadius"){                 scaleToRadius = Util::stob(argv[++i]);              continue; }
		if(arg=="--gradient" || arg=="-g"){         gradient = atoi(argv[++i]);                         continue; }
		if(arg=="--collisionFactor" || arg=="-c"){  collisionFactor = atof(argv[++i]);                  continue; }
		if(arg=="--decreaseSteps" ){                decreaseSteps = atoi(argv[++i]);                    continue; }
		if(arg=="--decreaseFactor"){                decreaseFactor = atof(argv[++i]);                   continue; }
		if(arg=="--stepSize"){                      stepSize = atof(argv[++i]);                         continue; }
		if(arg=="--maxRotation"   ){                maxRotation = atof(argv[++i]);                      continue; }
		if(arg=="--rejectsBeforeClose"){            poisson_max_rejects_before_close = atoi(argv[++i]); continue; }
		if(arg=="--metric"   ){                     metric_string = argv[++i];                          continue; }
		if(arg=="--planner"   ){                    planner_string = argv[++i];                         continue; }
		if(arg=="--rebuildLength"){                 rebuild_fragment_length = atoi(argv[++i]);          continue; }
		if(arg=="--rebuildFrequency"){              rebuild_frequency = atof(argv[++i]);                continue; }
		if(arg=="--rebuildInitial"){                rebuildInitialStructures = atoi(argv[++i]);         continue; }
		if(arg=="--rebuildAggression"){             rebuildAggression = atoi(argv[++i]);                continue; }
		if(arg=="--flexibleRibose"){                flexibleRibose = atoi(argv[++i]);                   continue; }
		if(arg=="--seed"){                          seed = atoi(argv[++i]);                             continue; }
		if(arg=="--biasToTarget" || arg=="-bias"){  biasToTarget = atof(argv[++i]);                     continue; }
    if(arg=="--convergeDistance"){              convergeDistance = atof(argv[++i]);                 continue; } //Previously rmsdThreshold
    if(arg=="--saveData"){                      saveData = atoi(argv[++i]);                         continue; }
		if(arg=="--sampleReverse"){                 sampleReverse = Util::stob(argv[++i]);              continue; }
    if(arg=="--residueNetwork" || arg=="-res"){ Selection::split( argv[++i],",",residueNetwork);    continue; }
    if(arg=="--alignAlways"){                   alignAlways = Util::stob(argv[++i]);                continue; } //Previously align
		if(arg=="--alignIni"){                      alignIni = Util::stob(argv[++i]);                   continue; }
		if(arg=="--selectAtoms"){                   selectAtoms = argv[++i];                            continue; }
		if(arg=="--preventClashes"){                preventClashes = Util::stob(argv[++i]);             continue; }
		if(arg=="--selectionAlign"){                selectionAlign = argv[++i];                         continue; }
		if(arg=="--selectionMoving"){               selectionMoving = argv[++i];                        continue; }
		if(arg=="--m_root"){                          root = atoi(argv[++i]);                             continue; }
    if(arg=="--projectConstraints"){            projectConstraints = Util::stob(argv[++i]);         continue; }
		if(arg=="--collisionCheck"){                collisionCheck = argv[++i];                         continue; }
		if(arg=="--frontSize"){                     frontSize = atoi(argv[++i]);                        continue; }
		if(arg=="--switchAfter"){                   switchAfter = atoi(argv[++i]);                      continue; }

		if(arg.at(0)=='-'){
			cerr<<"Unknown option: "<<arg<<endl<<endl;
			enableLogger("so");
			printUsage(argv[-1]);
			exit(-1);
		}
	}

	//Usage instructions if hbondFile is blank
	if( hydrogenbondMethod!="user" && hydrogenbondMethod!="dssr" && hydrogenbondMethod!="rnaview" && hydrogenbondMethod!="first" && hydrogenbondMethod!="FIRST" && hydrogenbondMethod!="vadar" ){
		enableLogger("so");
		cerr<<"Error: The hbond method ("<<hydrogenbondMethod<<") is not valid"<<endl;
		printUsage(argv[-1]);
		exit(-1);
	}
	if( hydrogenbondFile.empty() ){
		enableLogger("so");
		if(hydrogenbondMethod=="user"){
			log("so")<<"The hbond method '"<<hydrogenbondMethod<<"' was chosen.\nEach line in the file specified with --hbondFile should contain ";
			log("so")<<"two numbers indicating the atom-id of the hydrogen and the acceptor atoms in the initial structure file."<<endl;
		}else if(hydrogenbondMethod=="dssr"){
			log("so")<<"The hbond method '"<<hydrogenbondMethod<<"' was chosen.\nThe file specified with --hbondFile should contain the output from ";
			log("so")<<"the dssr program. There is a web-server where these files can be downloaded (http://w3dna.rutgers.edu/) ";
			log("so")<<"and dssr itself can be found at http://x3dna.org. To get the output (hbonds.txt) from a pdb-file (2B7G.pdb) ";
			log("so")<<"run the command\nx3dna-dssr -i=2B7G.pdb -o=hbonds.txt"<<endl;
			log("so")<<"Note that dssr only finds base-pair types, so exactly one hbond is added between U-A and G-C paired side-chains."<<endl;
		}else if(hydrogenbondMethod=="rnaview"){
			log("so")<<"The hbond method '"<<hydrogenbondMethod<<"' was chosen.\nThe file specified with --hbondFile should contain the output from ";
			log("so")<<"the RNAView program. The program can be downloaded from http://ndbserver.rutgers.edu/ndbmodule/services/download/rnaview.html ";
			log("so")<<"(possibly http://ndbserver.rutgers.edu/files/ftp/NDB/programs/rna/). To get the output (2B7G.pdb.out) from a pdb-file (2B7G.pdb) ";
			log("so")<<"run the command\nrnaview 2B7G.pdb"<<endl;
			log("so")<<"Note that rnaview only finds base-pair types, so exactly one hbond is added between U-A and G-C paired side-chains."<<endl;
		}else if(hydrogenbondMethod=="first" || hydrogenbondMethod=="FIRST"){
			log("so")<<"The hbond method '"<<hydrogenbondMethod<<"' was chosen.\nThe file specified with --hbondFile should contain the output from ";
			log("so")<<"the FIRST program generated by typing: FIRST <file>.pdb -non -hbout (the output will be named hbonds.out)"<<endl;
		}else if(hydrogenbondMethod=="vadar"){
			log("so")<<"The hbond method '"<<hydrogenbondMethod<<"' was chosen.\nThe file specified with --hbondFile should contain the output from ";
			log("so")<<"the vadar server (vadar.wishartlab.com) that has been edited so every line has this format (multiple spaces not important) ";
			log("so")<<"'92A     C      73A    N               1.91'"<<endl;
		}else{
			printUsage(argv[-1]);
			log("so")<<"No recognizable hbondMethod was specified ("<<hydrogenbondMethod<<")."<<endl;
			exit(-1);
		}
		exit(3);
	}else{
		if(hydrogenbondMethod.empty()) hydrogenbondMethod = "user"; //Default when a file is specified but no method
	}

	//Check initial structure
	if(initialStructureFile.empty()) {
		enableLogger("so");
		cerr<<"No initial structure file supplied"<<endl<<endl;
		exit(-1);
	}
	//Check target structure
	if(targetStructureFile.empty()){
		log("so")<<"No target structure file supplied, random sampling."<<endl<<endl;
	}

	//Check for valid metric
	std::transform(metric_string.begin(), metric_string.end(), metric_string.begin(), ::tolower);
	if(metric_string!="dihedral" && metric_string!="rmsd"){
		cerr<<"Invalid --m_metric option: "<<metric_string<<". Must be either 'dihedral' or 'rmsd'."<<endl;
		exit(-1);
	}

	//Check for valid planner
	std::transform(planner_string.begin(), planner_string.end(), planner_string.begin(), ::tolower);
	if( planner_string!="dihedralrrt" &&
      planner_string!="binnedrrt" &&
			planner_string!="dccrrt" &&
      planner_string!="poisson"){
		cerr<<"Invalid --planner option: "<<planner_string<<endl;
		exit(-1);
	}

  //Ensure that decreaseSteps>0 if preventClashes set
  if(preventClashes && decreaseSteps==0)
    decreaseSteps=5;

	//Set workingDirectory and moleculeName using the initialStructureFile.
	char* tmp = realpath(initialStructureFile.c_str(), nullptr);
	if(tmp==nullptr){ cerr<<initialStructureFile<<" is not a valid PDB-file"<<endl; exit(-1); }


	string pdb_file(tmp);
	int nameSplit = pdb_file.find_last_of("/\\");
	if(workingDirectory.empty()) {
    workingDirectory = pdb_file.substr(0, nameSplit + 1);
    log("so")<<"Changing working directory to "<<workingDirectory<<" using initial"<<endl;
  }

  //Check target structure
  if(targetStructureFile.empty()){
    log("so")<<"No target structure file supplied, random sampling."<<endl;
  }

  //Set default convergeDistance
  if(convergeDistance<0.0) {
    if(metric_string=="rmsd")     convergeDistance = 0.01;
    if(metric_string=="dihedral") convergeDistance = 1.0e-8;
  }

	//Search for hbonds file. Print warning if not found
	if(!fileExists(hydrogenbondFile)){
		string alt = workingDirectory+"/"+hydrogenbondFile;
		if(fileExists(alt))
			hydrogenbondFile = alt;
		else
			cerr<<"Warning: "<<hydrogenbondFile<<" not found. No hydrogen bonds added. Specify valid hydrogen bonds with the -hydrogens arg."<<endl;
	}

	//Error if rebuild is requested without annotations
	if( (annotationFile.empty() || !fileExists(annotationFile)) && hydrogenbondMethod!="rnaview" &&
	        (rebuild_frequency>0.0 || rebuildInitialStructures>0) ) {
		cerr<<"Error: No or invalid annotation file specified so the requested rebuilding can't be performed"<<endl;
		exit(-1);
	}

	//TODO: Exit with error if both FIRST and sugar conformations are enabled

	free(tmp);
}


void SamplingOptions::initializeVariables(){
	initialStructureFile    = "";
	targetStructureFile     = "";
	annotationFile          = "";
	hydrogenbondFile        = "";
	hydrogenbondMethod      = "";
	samplesToGenerate       = 10;
	explorationRadius       = 5.0;
  scaleToRadius           = false;
	sampleRandom            = true;
	gradient                = 0;
	collisionFactor         = 0.75;
	decreaseSteps 					= 0;
	decreaseFactor          = 0.5;
	stepSize                = 1.0;
	maxRotation             = 3.1415/18;
  poisson_max_rejects_before_close = 10;
	metric_string           = "rmsd";
	planner_string          = "binnedRRT";
	rebuild_fragment_length = 0;
	rebuild_frequency       = 0.0;
 rebuildInitialStructures = 0;
 rebuildAggression        = 0;
 flexibleRibose           = false;
 seed                     =  418;
 saveData                 =  0;
 sampleReverse            = false;
 biasToTarget             = 0.0;
 convergeDistance         = -1.0; ///<Changes depending on m_metric. Initialize <0, it is set depending on metric
 alignAlways              = false;
 alignIni                 = false;
 selectAtoms              = "heavy";
 preventClashes           = false;
 selectionAlign           = "";
 selectionMoving          = "";
 root                     = 0;
 projectConstraints       = true;
 collisionCheck           = "all";
 frontSize                = 50;
 switchAfter              = 20000;
}

void SamplingOptions::print(){
	log("so")<<"Sampling options:"<<std::endl;
	log("so")<<"\t--initial "<<initialStructureFile<<endl;
	log("so")<<"\t--target "<<targetStructureFile<<endl;
	log("so")<<"\t--annotation "<<annotationFile<<endl;
	log("so")<<"\t--hbondMethod "<<hydrogenbondMethod<<endl;
	log("so")<<"\t--hbondFile "<<hydrogenbondFile<<endl;
	log("so")<<"\t--extraCovBonds "; for(unsigned int i=0;i<extraCovBonds.size();i++) log("so")<<extraCovBonds[i]<<","; log("so")<<endl;
	log("so")<<"\t--workingDirectory "<<workingDirectory<<endl;
	//log("so")<<"\t--rigidbodiesFromFIRST "<<(rigidbodiesFromFIRST?"true":"false")<<endl;
	log("so")<<"\t--samples "<<samplesToGenerate<<endl;
	log("so")<<"\t--radius "<<explorationRadius<<endl;
	log("so")<<"\t--scaleToRadius "<<scaleToRadius<<endl;
	log("so")<<"\t--sampleRandom "<<sampleRandom<<endl;
	log("so")<<"\t--gradient "<<gradient<<endl;
	log("so")<<"\t--collisionFactor "<<collisionFactor<<endl;
	log("so")<<"\t--decreaseSteps "<< decreaseSteps <<endl;
	log("so")<<"\t--decreaseFactor "<<decreaseFactor<<endl;
	log("so")<<"\t--stepSize "<<stepSize<<endl;
	log("so")<<"\t--maxRotation "<<maxRotation<<endl;
  log("so")<<"\t--rejectsBeforeClose "<<poisson_max_rejects_before_close<<endl;
	log("so")<<"\t--metric "<<metric_string<<endl;
	log("so")<<"\t--planner "<<planner_string<<endl;
	log("so")<<"\t--rebuildLength "<<rebuild_fragment_length<<endl;
	log("so")<<"\t--rebuildFrequency "<<rebuild_frequency<<endl;
	log("so")<<"\t--rebuildInitial "<<rebuildInitialStructures<<endl;
	log("so")<<"\t--rebuildAggression "<<rebuildAggression<<endl;
	log("so")<<"\t--flexibleRibose "<<flexibleRibose<<endl;
	log("so")<<"\t--seed "<<seed<<endl;
	log("so")<<"\t--biasToTarget "<<biasToTarget<<endl;
  log("so")<<"\t--convergeDistance "<<convergeDistance<<endl;
	log("so")<<"\t--saveData "<<saveData<<endl;
	log("so")<<"\t--sampleReverse "<<sampleReverse<<endl;
	log("so")<<"\t--alignAlways "<<alignAlways<<endl;
	log("so")<<"\t--alignIni "<<alignIni<<endl;
	log("so")<<"\t--selectAtoms "<<selectAtoms<<endl;
	log("so")<<"\t--preventClashes "<<preventClashes<<endl;
	log("so")<<"\t--selectionAlign "<<selectionAlign<<endl;
	log("so")<<"\t--selectionMoving "<<selectionMoving<<endl;
	log("so")<<"\t--m_root "<<root<<endl;
	log("so")<<"\t--projectConstraints "<<projectConstraints<<endl;
	log("so")<<"\t--collisionCheck "<<collisionCheck<<endl;
	log("so")<<"\t--frontSize "<<frontSize<<endl;
	log("so")<<"\t--switchAfter "<<switchAfter<<endl;
}

void SamplingOptions::printUsage(char* pname){
	log("so")<<"Usage: "<<pname<<" explore [option list]"<<endl;
	log("so")<<"The KGS program will start sampling using the specified options."<<endl;
	log("so")<<"Options are:"<<endl;

	log("so")<<"\t--initial <path to structure> \t: Specifies the initial structure."<<endl;

	log("so")<<"\t--target <path to structure> \t: Specifies the target structure, optional."<<endl;

	log("so")<<"\t--annotation <file-path> \t: Annotations can specify secondary structures or other things ";
	log("so")<<"relevant to the sampling strategy. For RNA, standard WC will indicate non-free residues that wont be rebuilt"<<endl;

	log("so")<<"\t--hbondMethod <user|dssr|rnaview|first|vadar> \t: Format of the --hbondFile. If no --hbondFile argument is provided, instructions ";
	log("so")<<"how to generate a hbondFile are printed."<<endl;

	log("so")<<"\t--hbondFile <path to hydrogen bond file> \t: Hydrogen bond definition file. The format is specified by the choice ";
	log("so")<<"of --hbondMethod. Leave this field blank for instructions how to generate the hbond file."<<endl;
	log("so")<<"\t--extraCovBonds <resi1>/<name1>-<resi2>/<name2>[,...] \t: Extra covalent bonds. Can override an h-bond."<<endl;

	log("so")<<"\t--workingDirectory <directory> \t: Working directory. Output is stored here."<<endl;

	//log("so")<<"\t--rigidbodiesFromFIRST    \t: Indicates that FIRST should be used to identify rigid bodies. If not specified ";
	//log("so")<<"rigid bodies are identified using builtin definitions of partial double-bonds (works for RNA/DNA/proteins)"<<endl;

	log("so")<<"\t--samples, -s <whole number> \t: Indicates the number of samples to generate. Default is 10."<<endl;

	log("so")<<"\t--radius, -r <real number> \t: The exploration radius around the center structure. The default is 2Å."<<endl;

	log("so")<<"\t--sampleRandom  true/false>\t: if true, seed configurations are chosen at random. Default true."<<endl;

	log("so")<<"\t--scaleToRadius true/false>\t: if true, the random target configuration is scaled to inside the exploration radius. Default false."<<endl;

	log("so")<<"\t--gradient, -s <0|1|2|3|4> \t: Indicates method to calculate a new perturbation or gradient. 0 = random; 1 = torsion, no blending; 2 = torsion, low pass blending; 3 = msd, no blending; 4 = msd, low pass blending. Default is 0."<<endl;

	log("so")<<"\t--collisionFactor, -c <real number> \t: A number that is multiplied with the van der Waals radius when ";
	log("so")<<"checking for collisions. The default is 0.75."<<endl;

	//log("so")<<"\t--decreaseSteps <whole number> \t: If a non-colliding structure can not be found, try this many times to ";
	//log("so")<<"decrease the stepsize. Default is 1."<<endl;

	//log("so")<<"\t--decreaseFactor <real number> \t: If a non-colliding structure can not be found, decrease the stepsize ";
	//log("so")<<"by this factor. Default is 0.5."<<endl;

	log("so")<<"\t--stepSize <real number> \t: Initial step size to next sample as the norm of dihedral perturbation. Default 1."<<endl;

	log("so")<<"\t--maxRotation <real number> \t: The largest allowable change in torsion angle. The default is 10°."<<endl;

  log("so")<<"\t--rejectsBeforeClose <integer> \t: For poisson sampling: Number of perturbations attempted before closing. The default is 10°."<<endl;

	log("so")<<"\t--m_metric <rmsd|dihedral> \t: The m_metric to use in sampler. Default is rmsd."<<endl;

	log("so")<<"\t--planner <binnedRRT|dihedralRRT> \t: The planning strategy used to create samples. Default is binnedRRT."<<endl;

	log("so")<<"\t--rebuildLength <whole number>\t: The length of fragments that are rebuilt. Standard is 0."<<endl;

	log("so")<<"\t--rebuildFrequency <real number> \t: The frequency that rebuild perturbations will occur. The default is 0."<<endl;

	log("so")<<"\t--rebuildInitial <whole number>\t: The number of structures added to the initial pool by rebuilding free ";
	log("so")<<"loops. The default is 0."<<endl;

	log("so")<<"\t--rebuildAggression <0|1|2>\t: The aggression level of rebuilding. 0: Only sugars are resampled, backbone used for reclosing. ";

	log("so")<<"1: Sugars and side-chains are resampled, backbone used for reclosing. 2: Everything resampled, backbone used for reclosing."<<endl;

	log("so")<<"\t--flexibleRibose <0|1>\t: Indicate if ribose rings should be flexible (Standard: 1)"<<endl;

  log("so")<<"\t--seed <integer>\t: The seed used for random number generation (Standard: 418)"<<endl;

  log("so")<<"\t--biasToTarget, -bias <real number> \t: Percentage of using target as 'random seed configuration'. Default 0.1."<<endl;

  log("so")<<"\t--convergeDistance <real number> \t: The distance under which the goal conformation is considered reached. Default is 0.1 for RMSD and 1e-8 for Dihedral m_metric."<<endl;

  log("so")<<"\t--saveData <0|1|2>\t: Indicate whether files shall be saved! 0= none, 1=pdb and q, 2=all"<<endl;

  log("so")<<"\t--sampleReverse <true/false>\t: If true, iterative sampling from ini to goal and reverse"<<endl;

	log("so")<<"\t--alignIni "<<(alignIni?"true":"false: Align initial and target configuration in the beginning. Default false.")<<endl;

	log("so")<<"\t--selectAtoms <string atom selection> \t: Atom selection for gradients, distance computation etc., e.g. \"heavy\", \"name CA\", \"name C5\", \"backbone\", \"all\". Default heavy."<<endl;

	log("so")<<"\t--preventClashes "<<(preventClashes?"true":"false: Use clashing atoms to define additional constraints and prevent clash. Default true.")<<endl;

	log("so")<<"\t--selectionAlign <string of residue selection of molecule, e.g. \"resid 10 to 25 30 35 to 40\"> \t: Specifies the residues of the molecule that will undergo RMSD alignment during sampling."<<endl;

	log("so")<<"\t--selectionMoving <string of residue selection of molecule, e.g. \"resid 10 to 25 30 35 to 40\"> \t: Specifies the residues of the molecule that are used to determine the gradient."<<endl;

	log("so")<<"\t--m_root <integer>\t: The rigid body id for the m_root. Choose -1 to select the closest rigid body between m_molecule and target."<<endl;

	log("so")<<"\t--projectConstraints <true/false>\t: If false, then we don't project moves onto the constraint manifold. Only recommended for testing."<<endl;

	log("so")<<"\t--collisionCheck <string>\t: atoms used for collision detection: all (default), heavy, backbone"<<endl;

	log("so")<<"\t--frontSize <integer>\t: Size of the propagating front of samples in directed sampling."<<endl;

	log("so")<<"\t--switchAfter <integer>\t: Max number of steps before switching search directions (if bidirectional is active)."<<endl;

}



/** From http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c */
inline bool SamplingOptions::fileExists (const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

SamplingOptions* SamplingOptions::instance = nullptr;

SamplingOptions* SamplingOptions::getOptions()
{
  if(instance==nullptr) {
		cerr << "SamplingOptions::getInstance - Sampling options haven't been initialized"<<endl;
		exit(-1);
	}

  return instance;
}

SamplingOptions* SamplingOptions::createOptions(int argc, char* argv[] )
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new SamplingOptions(argc, argv);
  return instance;
}

SamplingOptions* SamplingOptions::createOptions()
{
  if(instance!=nullptr) {
    delete instance;
  }
  instance = new SamplingOptions();
  return instance;
}

//Todo: remove, but check in directions for usage
void SamplingOptions::setResidueNetwork(const Molecule * protein){//uses the selectionMoving and stores a "residueNetwork" list of ints
	//Determine if only certain residues shall be perturbed
	string selGradient = selectionMoving;
	Selection gradientSelection;
	vector<Atom*> atomsGradient;
	if(selGradient != ""){
		Selection gradientSelection(selGradient);
		vector<Residue*> residuesGradient = gradientSelection.getSelectedResidues(protein);
		for(vector<Residue*>::iterator ita = residuesGradient.begin(); ita != residuesGradient.end(); ++ita){
			residueNetwork.push_back( (*ita)->getId() );
//			log("dominik")<<"l "<<(*ita)->getId()<<endl;
		}
		log("dominik")<<" Number of residues for gradient = "<<residueNetwork.size()<<endl;
	}
}

void SamplingOptions::setAtomSets(const Molecule * protein, Molecule * target){
	//Here, we define the atom sets used to calculate a gradient, rmsd, or alignment.
	//If no input option is specified, then we use all atoms states in "atomsToChoose"
	//In case a target is present, only atoms are used from the selection that are present in both m_molecule structures

	//Within the user-provided selection of residues, we choose <selectAtoms> atoms (user-provided string, default "heavy")

	Selection alignSelection;
	vector<Atom*> atomsAlign;
	vector<Atom*> atomsMoving;

	if(selectionAlign != ""){
		Selection alignSelection(selectionAlign);
		vector<Residue*> residuesAlign = alignSelection.getSelectedResidues(protein);
		alignSelection.selection( selectAtoms ); // will use heavy atoms
		atomsAlign = alignSelection.getSelectedAtoms( residuesAlign );
	}
	else{
		Selection alignSelection(selectAtoms);
		atomsAlign = alignSelection.getSelectedAtoms(protein);
	}
	Atom *a1;
	std::string name, chainName;
	for (vector<Atom*>::iterator it=atomsAlign.begin(); it!=atomsAlign.end(); ++it) {
		a1=(*it);
		name = a1->getName();
		chainName = a1->getResidue()->getChain()->getName();
		int resId = a1->getResidue()->getId();
		if(target){
			Atom* a2=target->getAtom(chainName,resId, name);
			if(a2!=nullptr)
				m_atomsAlign.push_back(a1);
		} else{
			m_atomsAlign.push_back(a1);
		}
	}
	log("dominik")<<" Number atoms for alignment = "<<m_atomsAlign.size()<<endl;
//	for (vector<Atom*>::iterator it=m_atomsAlign.begin(); it!=m_atomsAlign.end(); ++it) {
//        log("dominik")<<"align res id: "<<(*it)->getResidue()->getId()<<", atom name: "<<(*it)->getName()<<endl;
//	}

	Selection movingSelection;

	if(selectionMoving != ""){
		Selection movingSelection(selectionMoving);
		vector<Residue*> residuesMoving = movingSelection.getSelectedResidues(protein);
		movingSelection.selection( selectAtoms ); // will use selected atoms only, default heavy
		atomsMoving = movingSelection.getSelectedAtoms( residuesMoving );
	}
	else{
		Selection movingSelection(selectAtoms);
		atomsMoving = movingSelection.getSelectedAtoms(protein);
	}
	for (vector<Atom*>::iterator it=atomsMoving.begin(); it!=atomsMoving.end(); ++it) {
		a1=(*it);
		name = a1->getName();
		chainName = a1->getResidue()->getChain()->getName();
		int resId = a1->getResidue()->getId();
		if(target){
			Atom* a2=target->getAtom(chainName,resId, name);
			if(a2!=nullptr)
				m_atomsMoving.push_back(a1);
		} else{
			m_atomsMoving.push_back(a1);
		}
	}
	log("dominik")<<" Number atoms for gradient = "<<m_atomsMoving.size()<<endl;
//	for (vector<Atom*>::iterator it=m_atomsMoving.begin(); it!=m_atomsMoving.end(); ++it) {
//        log("dominik")<<"gradient res id: "<<(*it)->getResidue()->getId()<<", atom name: "<<(*it)->getName()<<endl;
//	}
}
