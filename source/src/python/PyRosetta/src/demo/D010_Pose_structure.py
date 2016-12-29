#!usr/bin/env python

from __future__ import print_function

################################################################################
# A GENERAL EXPLANATION

"""
pose_structure.py

This script displays various structural and statistical contained within the
pose object. The Python syntax presented here is useful for quick investigation
structural data.

Instructions:

1) ensure that your PDB file is in the current directory
2) run the script:
    from commandline                        >python D010_Pose_structure.py

    from within python/ipython              [1]: run D010_Pose_structure.py

Author: Evan H. Baugh
    revised and motivated by Robert Schleif

Last updated by Boon Uranukul, 6/9/12

References:
    A. Leaver-Fay et al., "ROSETTA3: An object-oriented software suite for the
        simulation and design of macromolecules," Methods in Enymology 487,
        548-574 (2011).
"""

################################################################################
# THE BASIC PROTOCOL, pose_structure

"""
This sample script is setup for usage with
    commandline arguments,
    default running within a python interpreter,
    or for import within a python interpreter,
        (exposing all methods below)

The method pose_structure:
1.  obtains the pose's protein sequence
2.  obtains the pose's per-residue PDB number and icode
3.  obtains the pose's per-residue chain identification
4.  identifies the unique chain ids
5.  obtains the pose's secondary structure
6.  obtains the pose's per-residue backbone torsion angles
7.  outputs information on the requested residues

"""

import optparse    # for option sorting

from rosetta import *
from pyrosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

def pose_structure(pose, display_residues = []):
    """
    Extracts and displays various structural properties of the input  <pose>
        and its  <display_residues>  including:
            -PDB numbering
            -chain identification
            -sequence
            -secondary structure

    """
    # store the pose's number of residues, example Python syntax
    nres = pose.total_residue()

    # 1. obtain the pose's sequence
    sequence = pose.sequence()

    # 2. obtain a list of PDB numbering and icode as a single string
    pdb_info = pose.pdb_info()
    PDB_nums = [(str( pdb_info.number(i)) + pdb_info.icode(i)).strip()
        for i in range(1, nres + 1)]
    # 3. obtains a list of the chains organized by residue
    chains = [pdb_info.chain(i) for i in range(1, nres + 1)]
    # 4. extracts a list of the unique chain IDs
    unique_chains = []
    for c in chains:
        if c not in unique_chains:
            unique_chains.append(c)

    # start outputting information to screen
    print('\n' + '='*80)
    print('Loaded from' , pdb_info.name())
    print(nres , 'residues')
    print(len(unique_chains), 'chain(s) ('+ str(unique_chains)[1:-1] + ')')
    print('Sequence:\n' + sequence)

    # this object is contained in PyRosetta v2.0 and above
    # 5. obtain the pose's secondary structure as predicted by PyRosetta's
    #    built-in DSSP algorithm
    DSSP = protocols.moves.DsspMover()
    DSSP.apply(pose)    # populates the pose's Pose.secstruct
    ss = pose.secstruct()
    print( 'Secondary Structure:\n' + ss )
    print( '\t' + str(100. * ss.count('H') / len(ss))[:4] + '% Helical' )
    print( '\t' + str(100. * ss.count('E') / len(ss))[:4] + '% Sheet' )
    print( '\t' + str(100. * ss.count('L') / len(ss))[:4] + '% Loop' )

    # 6. obtain the phi, psi, and omega torsion angles
    phis = [pose.phi(i) for i in range(1, nres + 1)]
    psis = [pose.psi(i) for i in range(1, nres + 1)]
    omegas = [pose.omega(i) for i in range(1, nres + 1)]

    # this object is contained in PyRosetta v2.0 and above
    # create a PyMOLMover for exporting structures directly to PyMOL
    pymover = PyMOLMover()
    pymover.apply(pose)    # export the structure to PyMOL (optional)

    # 7. output information on the requested residues
    # use a simple dictionary to make output nicer
    ss_dict = {'L':'Loop', 'H':'Helix', 'E':'Strand'}
    for i in display_residues:
        print( '='*80 )
        print( 'Pose numbered Residue', i )
        print( 'PDB numbered Residue', PDB_nums[i-1] )
        print( 'Single Letter:', sequence[i-1] )
        print( 'Chain:', chains[i-1] )
        print( 'Secondary Structure:', ss_dict[ss[i-1]] )
        print( 'Phi:', phis[i-1] )
        print( 'Psi:', psis[i-1] )
        print( 'Omega:', omegas[i-1] )
        # extract the chis
        chis = [pose.chi(j + 1, i) for j in range(pose.residue(i).nchi() )]
        for chi_no in range(len(chis)):
            print( 'Chi ' + str(chi_no + 1) + ':', chis[chi_no] )
    print( '='*80 )

################################################################################
# INTERPRETING RESULTS

"""
This sample script is strictly to provide example syntax, it does not perform
any significant protocol and merely extracts data from a Pose. Please
investigate the accompanying scripts pose_scoring.py, fold_tree.py, movemap.py,
and packer_task.py to better understand the Pose object. The sample method and
PDB file presented for with this sample script work without error.
Since Rosetta is not entirely robust to all PDB files, several problems can
occur if the methods are modified or used with new PDB files. The most likely
problems could be:
    -if the input PDB file cannot be loaded into PyRosetta
    -if the residues to specifically investigate do not exist
"""

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls pose_structure with these values

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
# create an empty Pose object
pose = Pose()
# load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# default to the median residue number
residues = options.residues
if not options.residues:
    residues = [int(pose.total_residue()/2)]
elif options.residues == 'all':
    # accept the word 'all' in place of a residue list
    residues = range(1, pose.total_residue() + 1)
else:
    # please provide the residues of interest as, delimited
    residues = [int(r) for r in options.residues.split(',')]

pose_structure(pose, residues)

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
