
# Aligns model2 to model1 based on CA atoms and stores the result in `1cfc_model2_aligned.pdb`
kgs_align --initial 1cfc_model1.pdb --target 1cfc_model2.pdb --align "name CA"
