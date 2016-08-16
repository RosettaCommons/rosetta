#!/usr/bin/env python

################################################################################
# A GENERAL EXPLANATION

"""
fold_tree.py

This sample script explains the FoldTree object frequently changed in PyRosetta
protocols. It demonstrates how to set up a proper FoldTree for basic
applications and displays the structural impact of the FoldTree.

Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D030_Fold_tree.py

    from within python/ipython              [1]: run D030_Fold_tree.py

Author: Evan H. Baugh
    revised and motivated by Robert Schleif

Last updated by Boon Uranukul, 6/9/12

References:
    C. Wang, P. Bradley, and D. Baker, "Protein-protein docking with backbone
        flexibility," Jour. Mol. Bio. 373 (2) 503-519 (2007)

"""

################################################################################
# THE BASIC PROTOCOL, fold_tree

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing all methods below)

The method pose_parallel_objects:
1.  creates an example pose, simple short sequence
2.  setup values for the "middle" residues
3.  setup a FoldTree for modeling a single loop
        a.  using FoldTree.new_jump
        b.  using FoldTree.add_edge
4.  create a "linearized" version of the pose as an easy-to-see example,
        optionally write the coordinates of each change to a PDB file and
        export the structures to PyMOL
5.  successively alter the pose, display and (optionally) write the
        coordinates of each change to a PDB file
        a.  "linearize" the pose (linearized.pdb)
        b.  make a change before the first jump point (pre_jump.pdb)
        c.  make a change in the "first" jump edge (early_in_jump.pdb)
        d.  make a change in the "second" jump edge (late_in_jump.pdb)
        e.  make a change after the last jump point (post_jump.pdb)

"""

import optparse    # for option sorting

from rosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

def fold_tree(PDB_out = False):
    """
    Demonstrates the syntax necessary for basic usage of the FoldTree object
        performs these changes with a demonstrative pose example and writes
        structures to PDB files if  <PDB_out>  is True

    """

    ##########
    # FoldTree
    # a FoldTree encodes the internal coordinate dependencies of a Pose
    # a Pose object MUST have a FoldTree
    # the FoldTree allows regions of a pose to become independent,
    #    it is used in many important applications, particularly:
    #    loop modeling: where changes to the conformation of a loop region
    #    should NOT alter the conformation of the entire protein
    #    rigid-body docking: where changes in the position of one docking
    #    partner should not alter the position of the other docking partner
    # a FoldTree is effectively a list of Edge objects, you can view the Edges
    #    by printing the FoldTree ("print FoldTree")
    # the length of a FoldTree (FoldTree.size) MUST match the length of its
    #    corresponding Pose (Pose.total_residue)
    # it is possible to create an improper FoldTree, the method
    #    FoldTree.check_fold_tree returns True if the FoldTree is complete and
    #    usable and False otherwise
    # some Edge objects are Jumps, indicating a "jump" in the internal
    #    coordinate dependency
    # when a FoldTree is created, it can accept an optional integer argument
    #    setting the FoldTree to contain a single Edge with a length equal to
    #    the input integer value, the same result is attained by creating an
    #    empty FoldTree (no input arguments) and using the method
    #    FoldTree.simple_tree with an input integer equal to the size of the
    #    FoldTree

    # 1. create the example pose
    test_pose = pose_from_sequence('ACDEFGHIKLMNPQRSTVWY'*3)

    # 2. setup the jump points, where a jump is anchored, and the cutpoint
    cutpoint = test_pose.total_residue() / 2    # integer division, no decimal
    low_jump_point = cutpoint - 10
    high_jump_point = cutpoint + 10

    # the easiest way to create a new complete FoldTree is to use the method
    #    FoldTree.simple_tree to create and empty FoldTree and assign jumps to
    #    it using the method FoldTree.new_jump
    # the FoldTree constructor is overloaded to accept an input integer
    #    indicating how large to make the FoldTree

    # 3. create a simple, one jump FoldTree for the pose
    # a. using FoldTree.new_jump
    #pose_fold_tree = FoldTree(test_pose.total_residue())
    #### these two lines produce the same FoldTree as the one above
    pose_fold_tree = FoldTree()
    pose_fold_tree.simple_tree(test_pose.total_residue())

    pose_fold_tree.new_jump(low_jump_point, high_jump_point, cutpoint)
    print '\nThe first FoldTree is proper:', pose_fold_tree.check_fold_tree()

    # b. using FoldTree.add_edge
    # a more difficult method for creating a FoldTree is simply to create it
    #    empty and use the method FoldTree.add_edge to fill the FoldTree with
    #    new Edge data
    pose_fold_tree = FoldTree()
    pose_fold_tree.add_edge(1, low_jump_point, -1)
    pose_fold_tree.add_edge(low_jump_point, cutpoint, -1)
    pose_fold_tree.add_edge(low_jump_point, high_jump_point, 1)
    pose_fold_tree.add_edge(high_jump_point, test_pose.total_residue(), -1)
    pose_fold_tree.add_edge(high_jump_point, cutpoint + 1, -1)
    print 'The second FoldTree is proper:', pose_fold_tree.check_fold_tree()

    # demonstrate FoldTree's effect on structure
    # 4. linearize it
    for i in range(1, test_pose.total_residue() + 1):
        test_pose.set_phi(i, -180)
        test_pose.set_psi(i, 180)
        test_pose.set_omega(i, 180)

    # the Pose.fold_tree method is an overloaded getter/setter,
    #    providing it with no input returns the Pose's FoldTree object
    #    providing a FoldTree object as input overwrites the Pose's current
    #    FoldTree with the new one
    # the FoldTree is set here to prevent problems when "linearizing"
    test_pose.fold_tree(pose_fold_tree)

    # this object is contained in PyRosetta v2.0 and above (optional)
    pymover = PyMOL_Mover()

    # 5. change and display the new structures
    # a. export "linearized" structure
    test_pose.pdb_info().name('linearized')    # for PyMOL_Mover
    pymover.apply(test_pose)
    if PDB_out:
        test_pose.dump_pdb('linearized.pdb')
    print '\nlinearized structure output'

    # b. make an early change
    test_pose.set_phi(low_jump_point - 10, 50)
    test_pose.pdb_info().name('pre_jump')    # for PyMOL_Mover
    pymover.apply(test_pose)    # all downstream residues move
    if PDB_out:
        test_pose.dump_pdb('pre_jump.pdb')
    print 'pre jump perturbed structure output'

    # c. make a change in the first edge created by the jump
    test_pose.set_phi(low_jump_point + 5, 50)
    test_pose.pdb_info().name('early_in_jump')    # for PyMOL_Mover
    pymover.apply(test_pose)    # residues up to the cutpoint change
    if PDB_out:
        test_pose.dump_pdb('early_in_jump.pdb')
    print 'first internal jump edge perturbed structure output'

    # d. make a change in the second edge created by the jump
    test_pose.set_phi(high_jump_point - 5, 50)
    test_pose.pdb_info().name('late_in_jump')    # for PyMOL_Mover
    pymover.apply(test_pose)    # residues down to the cutpoint change
    if PDB_out:
        test_pose.dump_pdb('late_in_jump.pdb')
    print 'second internal jump edge perturbed structure output'

    # e. make a late change
    test_pose.set_phi(high_jump_point + 10, 50)
    test_pose.pdb_info().name('post_jump')    # for PyMOL_Mover
    pymover.apply(test_pose)    # all residues downstream move
    if PDB_out:
        test_pose.dump_pdb('post_jump.pdb')
    print 'post jump perturbed structure output'

################################################################################
# INTERPRETING RESULTS

"""
The PDB files output during the FoldTree section illustrate how the FoldTree
determines internal coordinate dependencies. This is easiest to observe these
changes by exporting these structures directly to PyMOL or by writing these
structures to PDB files (input PDB_out=True) and then loading these into PyMOL
or another viewer.
    load linearized.pdb
    load pre_jump.pdb
    load early_in_jump.pdb
    load late_in_jump.pdb
    load post_jump.pdb
If the PyMOL_Mover is utilized, the structures are stored in different PyMOL
objects. Cycle through them (starting with the connected straight
conformation) to observe the effects of a FoldTree.
The first is simple the "linearized" pose. Due to the logic of the method,
each change is applied successively so early_in_jump.pdb carries its
own change as well as that made for pre_jump.pdb (and so on). The PDB files
should thus be compared successively (ie. pre_jump.pdb to linearized.pdb ,
early_in_jump.pdb to pre_jump.pdb etc.)
Notice in pre_jump.pdb that every residue downstream (farther down the chain)
from the changed residue has moved as a result of the torsion change while
residues upstream remain unaffected.
Notice in early_in_jump.pdb that every residue downstream from the changed
residue and before the cutpoint has moved while all other residues remain
unaffected.
Notice in late_in_jump.pdb that every residue upstream from the changed residue
up to the cutpoint has moved while all other residues remain unaffected.
Notice in post_jump.pdb that every residue downstream from the changed residue
has moved while all other residues remain unaffected.
If the original "linearized" chain is not continuous (all adjacent residues
connected with a bond) the result may be more confusing.

"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls fold_tree with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb" with reduced
#    cycles/jobs to provide results quickly
parser = optparse.OptionParser()
parser.add_option('--PDB_out', dest = 'PDB_out',
    default = '',    # default to False
    help = 'flag for dumping structures to PDB files')
(options,args) = parser.parse_args()

# PDB_out flag
PDB_out = bool(options.PDB_out)

fold_tree(PDB_out)
