#!/bin/bash

ADDMTSSL=../../scripts/add_mtssl.py

# First argument is pdb-file, second is residue number to change
python $ADDMTSSL 1crn.pdb 23 > 1crn_mtssl.pdb
