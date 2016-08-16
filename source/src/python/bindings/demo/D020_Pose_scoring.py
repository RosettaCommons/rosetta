#!/usr/bin/env python

################################################################################
# A GENERAL EXPLANATION

"""
pose_scoring.py

This script displays methods for scoring pose objects and extracting score date.
The Python syntax presented here is useful for quick investigation of
Rosetta scores, including numerous statistical properties.

Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D020_Pose_scoring.py

    from within python/ipython              [1]: run D020_Pose_scoring.py

Author: Evan H. Baugh
    revised and motivated by Robert Schleif

Last updated by Boon Uranukul, 6/9/12

References:
    A. Leaver-Fay et al., "ROSETTA3: An object-oriented software suite for the
        simulation and design of macromolecules," Methods in Enymology 487,
        548-574 (2011).

"""

################################################################################
# THE BASIC PROTOCOL, pose_scoring

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing all methods below)

The method pose_scoring:
1.  scores the pose using the default fullatom (fa) ScoreFunction by:
        a.  using get_fa_scorefxn
        b.  using create_score_function
        c.  using create_score_function_ws_patch
        d.  using create_score_function and patching it after-wards
        e.  creating an empty ScoreFunction and manually setting the weights
2.  outputs the ScoreFunction evaluations
3.  obtains the pose's per-residue total scores
4.  obtains the non-zero weights of the ScoreFunction, the active ScoreTypes
5.  obtains the pose's per-residue energies for all active ScoreTypes
6.  obtains the hydrogen bond information
7.  creates a dictionary for per-residue hydrogen bond information
8.  approximates the pose's radius of gyration
9.  outputs the pose information
10. outputs information on the requested residues

"""

import optparse    # for option sorting

from rosetta import *
from rosetta.PyMolLink import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

def pose_scoring(pose, display_residues = []):
    """
    Extracts and displays various score evaluations on the input  <pose>
		and its  <display_residues>  including:
		total score
		score term evaluations
		per-residue score
		radius of gyration (approximate)

    """
    # this object is contained in PyRosetta v2.0 and above
    # create a PyMOL_Mover for exporting structures directly to PyMOL
    pymover = PyMOL_Mover()

    # 1. score the pose using the default full-atom (fa) ScoreFunction
    # a ScoreFunction is essentially a list of weights indicating which
    #    ScoreTypes are "on" and how they are scaled when summed together
    # ScoreTypes define what scoring methods are performed on the pose
    # some ScoreTypes are centroid or fullatom specific, others are not
    # there are (at least) 5 ways to obtain the general fullatom ScoreFunction
    # the 1st (a) and 2nd (b) methods are only present since PyRosetta v2.0
    #
    # a. this method returns the hard-coded set of weights for the standard
    #    fullatom ScoreFunction, which is currently called "talaris2013"
    fa_scorefxn = get_fa_scorefxn()

    # b. this method returns a ScoreFunction with its weights set based on
    #    files stored in the database/scoring/weights (.wts files)
    full_scorefxn = create_score_function('talaris2013')

    # c. this method sets the weights based on 'talaris2013.wts' and then
    #    corrects them based on 'docking.wts_patch'
    ws_patch_scorefxn = create_score_function_ws_patch('talaris2013', 'docking')

    # d. this method returns a ScoreFunction with its weights set by loading
    #    weights from 'talaris2013' followed by an adjustment by setting
    #    weights from 'docking.wts_patch'
    patch_scorefxn = create_score_function('talaris2013')
    patch_scorefxn.apply_patch_from_file('docking')

    # e. here an empty ScoreFunction is created and the weights are set manually
    scorefxn = ScoreFunction()
    scorefxn.set_weight(fa_atr, 0.800)    # full-atom attractive score
    scorefxn.set_weight(fa_rep, 0.440)    # full-atom repulsive score
    scorefxn.set_weight(fa_sol, 0.750)    # full-atom solvation score
    scorefxn.set_weight(fa_intra_rep, 0.004)    # f.a. intraresidue rep. score
    scorefxn.set_weight(fa_elec, 0.700)    # full-atom electronic score
    scorefxn.set_weight(pro_close, 1.000)    # proline closure
    scorefxn.set_weight(hbond_sr_bb, 1.170)    # short-range hbonding
    scorefxn.set_weight(hbond_lr_bb, 1.170)    # long-range hbonding
    scorefxn.set_weight(hbond_bb_sc, 1.170)    # backbone-sidechain hbonding
    scorefxn.set_weight(hbond_sc, 1.100)    # sidechain-sidechain hbonding
    scorefxn.set_weight(dslf_fa13, 1.000)    # disulfide full-atom score
    scorefxn.set_weight(rama, 0.200)    # ramachandran score
    scorefxn.set_weight(omega, 0.500)    # omega torsion score
    scorefxn.set_weight(fa_dun, 0.560)    # fullatom Dunbrack rotamer score
    scorefxn.set_weight(p_aa_pp, 0.320)
    scorefxn.set_weight(ref, 1.000)    # reference identity score

    # ScoreFunction a, b, and e above have the same weights and thus return
    #    the same score for an input pose.  Likewise, c and d should return the
    #    same scores.
    # 2. output the ScoreFunction evaluations
    #ws_patch_scorefxn(pose)    # to prevent verbose output on the next line
    print '='*80
    print 'ScoreFunction a:', fa_scorefxn(pose)
    print 'ScoreFunction b:', full_scorefxn(pose)
    print 'ScoreFunction c:', ws_patch_scorefxn(pose)
    print 'ScoreFunction d:', patch_scorefxn(pose)
    print 'ScoreFunction e:', scorefxn(pose)
    pose_score = scorefxn(pose)

    # 3. obtain the pose Energies object and all the residue total scores
    energies = pose.energies()
    residue_energies = [energies.residue_total_energy(i)
        for i in range(1, pose.total_residue() + 1)]

    # 4. obtain the non-zero weights of the ScoreFunction, active ScoreTypes
    weights = [ScoreType(s)
        for s in range(1, end_of_score_type_enumeration + 1)
        if scorefxn.weights()[ScoreType(s)]]
    # 5. obtain all the pose energies using the weights list
    # Energies.residue_total_energies returns an EMapVector of the unweighted
    #    score values, here they are multiplied by their weights
    #    remember when performing individual investigation, these are the raw
    #    unweighted score!
    residue_weighted_energies_matrix = [
        [energies.residue_total_energies(i)[w] * scorefxn.weights()[w]
        for i in range(1, pose.total_residue() + 1)]
        for w in weights]

    # Unfortunately, hydrogen bonding scores are NOT stored in the structure
    #    returned by Energies.residue_total_energies
    # 6. hydrogen bonding information must be extracted separately
    pose_hbonds = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set( pose , False , pose_hbonds )

    # 7. create a dictionary with the pose residue numbers as keys and the
    #    residue hydrogen bonding information as values
    # hydrogen bonding information is stored as test in the form:
    # (donor residue) (donor atom) => (acceptor residue) (accecptor atom) |score
    hbond_dictionary = {}
    for residue in range(1, pose.total_residue() + 1):
        hbond_text = ''
        for hbond in range(1, pose_hbonds.nhbonds() + 1):
            hbond = pose_hbonds.hbond(hbond)
            acceptor_residue = hbond.acc_res()
            donor_residue = hbond.don_res()
            if residue == acceptor_residue or residue == donor_residue:
                hbond_text += str(donor_residue).ljust(4) + ' ' + \
                    str(pose.residue(donor_residue).atom_name(\
                        hbond.don_hatm() )).strip().ljust(4) + \
                    ' => ' + str(acceptor_residue).ljust(4) + ' ' + \
                    str(pose.residue(acceptor_residue).atom_name(\
                        hbond.acc_atm() )).strip().ljust(4) + \
                    ' |score: ' + str(hbond.energy()) + '\n'
        hbond_dictionary[residue] = hbond_text

    # 8. approximate the radius of gyration
    # there is an rg ScoreType in PyRosetta for performing this computation so
    #    a ScoreFunction can be made as an Rg calculator, likewise you can
    #    bias a structure towards more or less compact structures using this
    # NOTE: this is NOT the true radius of gyration for a protein, it uses
    #    the Residue.nbr_atom coordinates to save time, this nbr_atom is the
    #    Residue atom closest to the Residue's center of geometry
    RadG = ScoreFunction()
    RadG.set_weight(rg , 1)
    pose_radg = RadG(pose)

    # 9. output the pose information
	# the information is not expressed sequentially as it is produced because
	#    several PyRosetta objects and methods output intermediate information
	#    to screen, this would produce and unattractive output
    print '='*80
    print 'Loaded from' , pose.pdb_info().name()
    print pose.total_residue() , 'residues'
    print 'Radius of Gyration ~' , pose_radg
    print 'Total Rosetta Score:' , pose_score
    scorefxn.show(pose)
    # this object is contained in PyRosetta v2.0 and above
    pymover.apply(pose)
    pymover.send_energy(pose)

    # 10. output information on the requested residues
    for i in display_residues:
        print '='*80
        print 'Pose numbered Residue' , i
        print 'Total Residue Score:' , residue_energies[i-1]
        print 'Score Breakdown:\n' + '-'*45
        # loop over the weights, extract the scores from the matrix
        for w in range(len(weights)):
            print '\t' + name_from_score_type(weights[w]).ljust(20) + ':\t' ,\
                residue_weighted_energies_matrix[w][i-1]
        print '-'*45
        # print the hydrogen bond information
        print 'Hydrogen bonds involving Residue ' + str(i) + ':'
        print hbond_dictionary[i][:-1]
    print '='*80

################################################################################
# INTERPRETING RESULTS

"""
This sample script is strictly to provide example syntax, it does not perform
any significant protocol and merely extracts data from a Pose and from scoring
information. The sample method and PDB file presented with this sample script
work without error. Since Rosetta is not entirely robust to all PDB files,
several problems can occur if the methods are modified or used with new
PDB files. The most likely problems could be:
    -if the input PDB file cannot be loaded into PyRosetta
    -if the residues to specifically investigate do not exist
    -if the input pose is not fullatom

"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls pose_scoring with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb" with reduced
#    cycles/jobs to provide results quickly
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '../test/data/test_in.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel')
parser.add_option('--residues', dest = 'residues',
    default = '',    # default to the median residue number
    help = 'the (pose numbered) residues to inspect carefully')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
# create a pose from the desired PDB file
# -create an empty Pose object
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# default to the median residue number
residues = options.residues
if not options.residues:
    residues = [pose.total_residue()/2]
elif options.residues == 'all':
    # accept the word 'all' in place of a residue list
    residues = range(1, pose.total_residue() + 1)
else:
    # please provide the residues of interest as , delimited
    residues = [int(r) for r in options.residues.split(',')]

pose_scoring(pose, residues)

################################################################################
# ALTERNATE SCENARIOS

#################################
# Obtaining and Editing PDB files
"""
PDB files are the keys to structural Bioinformatics and structure prediction.
PDB files are most easily obtained from the RCSB website but may contain
variability which makes them incompatible with PyRosetta. To obtain a new
PDB file:

1) locate your protein of interest at http://www.pdb.prg/
2) download the PDB file, using a browser this includes:
    a. clicking "Download Files" on the upper right
    b. clicking "PDB File (text)", the second option
3) Manually edit the file to remove lines which may hinder PyRosetta
    (use PyMOL, grep, awk, Python, Biopython, or whatever technique you prefer)

Methods for downloading and generically "cleaning" PDB files should accompany
future PyRosetta releases.

"""
