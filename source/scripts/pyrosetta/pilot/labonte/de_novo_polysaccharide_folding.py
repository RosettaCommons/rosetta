#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
"""Brief:   This PyRosetta script folds a polysaccharide chain.

Author:  Jason W. Labonte

"""
# Imports
from os import mkdir
from argparse import ArgumentParser

from rosetta import init, get_fa_scorefxn, create_score_function, \
                    standard_packer_task, Pose, MoveMap, RotamerTrialsMover, \
                    SmallMover, ShearMover, MonteCarlo, PyMOL_Mover, \
                    hbond_sr_bb
from rosetta.core.pose import pose_from_saccharide_sequence
from structural_polysaccharide_refinement import refine_saccharide


# Parse arguments.
parser = ArgumentParser(description=__doc__)
parser.add_argument('n_residues', type=int,
                    help='the number of residues')
parser.add_argument('monosaccharide_name',
                    help='the abbreviated name of a single monosaccharide' +
                    ' unit, e.g., beta-D-Glcp')
parser.add_argument('linkage', type=str, default='(1->4)',
                    help='the linkage designation, in quotes, with the' +
                    ' following format: "(1->4)"')
parser.add_argument('--mute', action='store_true',
					help='flag to mute output for cycles')
parser.add_argument('--make_movie', action='store_true',
                    help='flag to output a movie frame for each step')
parser.add_argument('--pm', action='store_true',
                    help='flag to observe the protocol in PyMOL')
parser.add_argument('--mm', action='store_true',
                    help='flag to use the molecular mechanics score function')
parser.add_argument('--n_cycles', type=int, default=100,
                    help='the number of Monte Carlo cycles')
parser.add_argument('--kt', type=float, default=1.0,
                    help='the "temperature" to use for the MC cycles')
parser.add_argument('--n_small_moves', type=int, default=5,
                    help='the number of small moves for each MC cycle')
parser.add_argument('--n_shear_moves', type=int, default=5,
                    help='the number of shear moves for each MC cycle')
parser.add_argument('--max_angle', type=float, default=30.0,
                    help='the maximum angle for a backbone move')
parser.add_argument('--Hbond_weight', type=float, default=2.0,
                    help='the scoring weight for hydrogen bonds')
parser.add_argument('--min_type', default='linmin',
                    choices=['linmin', 'dfpmin'],
                    help='the type of minimization')
args = parser.parse_args()

# Initialize Rosetta.
init(extra_options='-include_sugars -override_rsd_type_limit -mute all')

# Format sequence.
sequence = args.monosaccharide_name + ("-" + args.linkage + "-" + \
                                args.monosaccharide_name) * (args.n_residues-1)

# Create pose.
print 'Creating pose with sequence:', sequence
pose = pose_from_saccharide_sequence(sequence)

# Set filenames.
output_base_filename = args.monosaccharide_name + '-' + \
                                                  str(args.n_residues) + '-mer'

# Run protocol.
refine_saccharide(pose, args, output_base_filename)