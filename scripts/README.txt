This folder stores python scripts used for post-processing. You can individually rerun scripts from the command line via python ../scriptToRun.py <inputs> and adapt input parameters as desired. The current script list provides an extensive number of graphs and output to create publication-ready figures in pyMol.

Current post-processing that will automatically be generated upon simulation:

- RMSD of all samples to initial and target (plotBasics.py)
- pairwise-RMSD plots for the bidirectional exploration trees 
- vdW energy along the paths

- constraint-related plots: bond differences between initial and target, h-bond violations along the paths (plotConstraints.py, plotBondLengths.py)
- multi-model pdb files to create „movies“ in pymol (default step: 10, combineKGSpath_steps.py)
- insights into vdW energy and clashes (pathPlots.py)
- pdb / pml files for torsion RMSF cartoon putty representations in pyMol (computeAbsoluteSteps.py, colorByTorsionRMSF.py)
- fancy forward and reverse exploration tree with encountered collisions in color (pairwiseTreeDistances.py)
- three different representations for clash analysis: clash density distribution on the molecule, residue network graphs, and disjoint clash networks for pymol (all clash… .py files)

You can find more details in the individual Python Scripts and the ../experiments/Makefile_inc file.

Happy sampling!