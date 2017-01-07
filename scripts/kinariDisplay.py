#!/usr/bin/env python

import __main__
__main__.pymol_argv = [ 'pymol', '-q'] # Quiet and no GUI

import sys, time, os
try:
    import pymol
except:
    print "The pymol site-packages must be in your PYTHONPATH"
    sys.exit(-1)

pymol.finish_launching()

structure_name = None

def read_bonds(fname):
    print "read_bonds(",fname,")"
    ret = []
    with open(fname) as fhandle:
        for line in fhandle:
            print line
            tokens = line.split()
            if len(tokens)<5: continue
            ret.append( (tokens[1], tokens[3], float(tokens[4])) )
    return ret

hbonds = []
rbonds = []
ccbonds = []

for arg in sys.argv:
    if arg.endswith(".pdb") or arg.endswith(".pdb.knr"):
        if structure_name: 
            print "Trying to read",arg,", but a structure was already read:",structure_name+", Please only read 1 PDB"
            pymol.cmd.quit()
            sys.exit(-1)

        structure_name = arg.split('/')[-1].split('.')[0]
        pymol.cmd.load(arg, structure_name)

    if arg.lower().endswith("hbonds.bnd.knr") or arg.lower().endswith("hbonds.pruned.bnd.knr"):
        hbonds = read_bonds(arg)

    if arg.lower().endswith("rbonds.bnd.knr") or arg.lower().endswith("rbonds.pruned.bnd.knr"):
        rbonds = read_bonds(arg)

    if arg.lower().endswith("disulfidebonds.bnd.knr"):
        ccbonds = read_bonds(arg)

pymol.cmd.set_color("slate",[0.7,0.8,0.9])

pymol.cmd.set("dash_gap", "0")

print "HBONDS:",hbonds
# Display hydrogen bonds
if hbonds:
    for (id1,id2,energy) in hbonds:
        pymol.cmd.distance("hbonds", "(ID "+id1+")", "(ID "+id2+")")

    pymol.cmd.show("everything", "hbonds")
    pymol.cmd.hide("labels", "hbonds")
    pymol.cmd.color("slate", "hbonds")

# Display resonance bonds
if rbonds:
    for (id1,id2,energy) in rbonds:
        pymol.cmd.distance("rbonds", "(ID "+id1+")", "(ID "+id2+")")

    pymol.cmd.show("everything", "rbonds")
    pymol.cmd.hide("labels", "rbonds")
    pymol.cmd.color("red", "rbonds")
    pymol.cmd.disable("rbonds")



# Display disulphide bonds
if ccbonds:
    for (id1,id2,energy) in ccbonds:
        pymol.cmd.distance("disulphide-bonds", "(ID "+id1+")", "(ID "+id2+")")

    pymol.cmd.show("everything", "disulphide-bonds")
    pymol.cmd.hide("labels", "disulphide-bonds")
    pymol.cmd.color("yellow", "disulphide-bonds")






