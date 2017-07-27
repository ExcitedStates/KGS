"""
Compares coordinates of KGS ensemble to MD ensemble
Calculates average RMSD between each KGS member and closest RMSD member (a)
Calculates maximum RMSD within MD ensemble (b)
Returns fraction of a/b
Strong overlap if less than one, dissimilar if greater than one
"""


from pdb_structure import *
import itertools


def getMinRMSD(kgs_model, md_ensemble):
    """Returns RMSD between given kgs member and closest MD member."""
    min_rmsd = kgs_model.rmsd(md_ensemble.models[0], ['CA'])
    for md_model in md_ensemble.models:
        rmsd = kgs_model.rmsd(md_model, ['CA'])  # If source of error: use pdb files, and reference models
        if rmsd < min_rmsd:
            min_rmsd = rmsd
    return min_rmsd


def getAvgMinRMSD(kgs_ensemble, md_ensemble):
    """Returns average RMSD between each kgs member and closest md member."""
    rmsd_vals = []
    for kgs_member in kgs_ensemble.models:
        rmsd = getMinRMSD(kgs_member, md_ensemble)
        rmsd_vals.append(rmsd)
    avg_rmsd = sum(rmsd_vals)/len(rmsd_vals)
    return avg_rmsd


def getMaxRMSD(md_ensemble):
    """Returns maximum RMSD in given ensemble."""
    max_rmsd = 0
    for x, y in itertools.combinations(md_ensemble.models, 2):
        rmsd = x.rmsd(y, ['CA'])
        if rmsd > max_rmsd:
            max_rmsd = rmsd
    return max_rmsd


kgs_ensemble = PDBFile("KGS_Results.pdb")
md_ensemble = PDBFile("MD_Results.pdb")
min_kgs = getAvgMinRMSD(kgs_ensemble, md_ensemble)
max_md = getMaxRMSD(md_ensemble)
print("Score : " + str(max_md/min_kgs))
