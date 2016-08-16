#!/usr/bin/env python

################################################################################
# A GENERAL EXPLANATION

"""
movemap.py

This sample script explains the MoveMap object frequently changed in PyRosetta
protocols. It demonstrates how to set up a proper MoveMap for basic
applications and applies this on backbone torsion minimization.

Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D040_Movemap.py

    from within python/ipython              [1]: run D040_Movemap.py

Author: Evan H. Baugh
    revised and motivated by Robert Schleif

Last updated by Boon Uranukul, 6/9/12

References:
    C. Wang, P. Bradley, and D. Baker, "Protein-protein docking with backbone
        flexibility," Jour. Mol. Bio. 373 (2) 503-519 (2007)

"""

# WARNING
"""
Yes, this script should (conventionally) be called move_map.py however the
MoveMap object is frequently implemented against this convention (for example
MinMover.movemap) and as such, I decided to continue this tradition since it
matches the syntax.

"""

################################################################################
# THE BASIC PROTOCOL, movemap

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing all methods below)

The method movemap:
1.  setup a MoveMap for minimizing the pose score by altering the backbone
        torsions of a set of residues
2.  perform minimization
        -display the change in score
        -write the minimized pose to a PDB file (minimized.pdb)

"""

import optparse    # for option sorting

from rosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


def movemap(pose, PDB_out = False):
    """
    Demonstrates the syntax necessary for basic usage of the MoveMap object
        performs these changes with a demonstrative backbone minimization
        using  <pose>  and writes structures to PDB files if  <PDB_out>  is True

    """

    #########
    # MoveMap
    # a MoveMap encodes what data is allowed to change in a Pose, referred to as
    #    its degrees of freedom
    # a MoveMap is separate from a Pose and is usually required by a Mover so
    #    that the correct degrees of freedom are manipulated, in this way,
    #    MoveMap and Pose objects often work in parallel
    # several MoveMap's can correspond to the same Pose
    # a MoveMap stores information on a per-residue basis about the
    #    backbone ({phi, psi, omega}) and chi ({chi_i}) torsion angle sets
    #    the MoveMap can only set these sets of torsions to True or False, it
    #    cannot set freedom for the individual angles (such as phi free and psi
    #    fixed)
    # the MoveMap has no upper-limit on its residue information, it defaults to
    #    all residues (up to residue 99999999) backbone and chi False
    # you can view the MoveMap per-residue torsion settings by using the
    #    MoveMap.show( Pose.total_residue() ) method (the input argument is the
    #    highest residue to output, it does not support viewing a range)
    pose_move_map = MoveMap()
    # change all backbone torsion angles
    pose_move_map.set_bb(True)
    # change all chi angle torsion angles (False by default)
    pose_move_map.set_chi(False)
    # change a single backbone torsion angles
    #pose_move_map.set_bb(1, True)    # example syntax
    # change a single residue's chi torsion angles
    #pose_move_map.set_chi(1, True)    # example syntax
    pose_move_map.show(pose.total_residue())

    # perform gradient based minimization on the "median" residues, this
    #    method (MinMover) determines the gradient of an input pose using a
    #    ScoreFunction for evaluation and a MoveMap to define the degrees of
    #    freedom
    # create a standard ScoreFunction
    scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')
    # redefine the MoveMap to include the median half of the residues
    # turn "off" all backbone torsion angles
    pose_move_map.set_bb(False)    # reset to backbone False
    # turn "on" a range of residue backbone torsion angles
    pose_move_map.set_bb_true_range(pose.total_residue() / 4,
        pose.total_residue() * 3 / 4)
    # create the MinMover
    minmover = MinMover()
    minmover.score_function(scorefxn)
    minmover.movemap(pose_move_map)

    # create a copy of the pose
    test_pose = Pose()
    test_pose.assign(pose)
    # apply minimization
    scorefxn(test_pose)    # to prevent verbose output on the next line

    pymover = PyMOL_Mover()
    #### uncomment the line below and "comment-out" the two lines below to
    ####    export the structures into different PyMOL states of the same object
    #pymover.keep_history = True    # enables viewing across states

    #### comment-out the line below, changing PDBInfo names tells the
    ####    PyMOL_Mover to produce new objects
    test_pose.pdb_info().name('original')
    pymover.apply(test_pose)
    print '\nPre minimization score:', scorefxn(test_pose)

    minmover.apply(test_pose)
    if PDB_out:
        test_pose.dump_pdb('minimized.pdb')

    print 'Post minimization score:', scorefxn(test_pose)
    #### comment-out the line below
    test_pose.pdb_info().name('minimized')
    pymover.apply(test_pose)

################################################################################
# INTERPRETING RESULTS

"""
Cycle through the structures exported to PyMOL to observe the structural change
performed by backbone minimization.
The PDB file output demonstrating minimization and MoveMap degree of freedom
selection should be compared to the original PDB file (test_in.pdb for the
sample case).
Notice that the backbone residues of the "middle" section have moved slightly.
Since only the backbone changed, PyMOL's cartoon representation (and other
similar views) displays this change well. Even though only the "middle" residue
torsions have changed, the downstream residue coordinates have moved while the
upstream residue coordinates remain stationary. Try comparing the minimization
using successive states of the same PyMOL object or by loading in the PDB files
separately and aligning them.

"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls movemap with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb"
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '../test/data/test_in.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel')
parser.add_option('--PDB_out', dest = 'PDB_out',
    default = '',    # default to False
    help = 'flag for dumping structures to PDB files')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
# create a pose from the desired PDB file
# -create an empty Pose object
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# PDB_out flag
PDB_out = bool(options.PDB_out)

movemap(pose, PDB_out)
