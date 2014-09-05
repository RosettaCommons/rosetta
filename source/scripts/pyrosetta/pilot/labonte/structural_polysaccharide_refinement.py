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
                    SmallMover, ShearMover, MinMover, MonteCarlo, \
                    PyMOL_Mover, hbond_sr_bb, hbond_bb_sc, hbond_sc, \
                    fa_elec, pose_from_pdb
from rosetta.core.chemical.carbohydrates import CarbohydrateInfo
from rosetta.protocols.simple_moves import RingConformationMover


def refine_saccharide(pose, args, output_base_filename):
    """Perform packing, small, and shear moves."""
    # Set up ScoreFunction.
    if args.mm:
        sf = create_score_function('mm_std')
    else:
        sf = get_fa_scorefxn()
    sf.set_weight(hbond_sr_bb, args.Hbond_weight)
    sf.set_weight(hbond_bb_sc, args.Hbond_weight)
    sf.set_weight(hbond_sc, args.Hbond_weight)
    sf.set_weight(fa_elec, 2.0)

    print ' Initial Score:', sf(pose)
    sf.show(pose)  # TEMP

    # Set up packer.
    pt = standard_packer_task(pose)
    pt.restrict_to_repacking()
    pt.or_include_current(True)

    packer = RotamerTrialsMover(sf, pt)

    # Set up BB Movers.
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_nu(False)

    small_mover = SmallMover(mm, args.kt, args.n_small_moves)
    small_mover.angle_max(args.max_angle)
    shear_mover = ShearMover(mm, args.kt, args.n_shear_moves)
    shear_mover.angle_max(args.max_angle/2)
    min_mover = MinMover(mm, sf, args.min_type, 0.001, True)

    # Set up ring Mover.
    if not args.lock_rings:
        mm_rings = MoveMap()
        mm_rings.set_bb(False)
        mm_rings.set_chi(False)
        mm_rings.set_nu(True)

        ring_flipper = RingConformationMover()
        ring_flipper.movemap(mm_rings)

    # Set up Monte Carlo object.
    mc = MonteCarlo(pose, sf, args.kt)

    # Set filenames.
    output_pdb_filename = output_base_filename + '.pdb'

    if args.make_movie:
        output_frame_base_filename = output_base_filename + '_frame'
        dir_name = output_base_filename + '_movie_frames'
        try:
            mkdir(dir_name)
        except OSError:
            print 'Warning:',
            print 'Directory already exists; files will be overwritten.'

    if args.pm:
        pm = PyMOL_Mover()

    # Prepack.
    packer.apply(pose)
    if args.pm:
        pm.apply(pose)
    print ' Score after prepacking:', sf(pose)
    sf.show(pose)  # TEMP

    # Fold the polysaccharide.
    print 'Folding/Refining structural polysaccharide with',
    print pose.total_residue(),
    #print pose.residue(2).carbohydrate_info().short_name(), 'residues...'
    print 'residues...'

    for cycle in range(1, args.n_cycles + 1):
        small_mover.apply(pose)

        shear_mover.apply(pose)

        if not args.lock_rings:
            ring_flipper.apply(pose)

        packer.apply(pose)

        min_mover.apply(pose)

        mc.boltzmann(pose)
        if not args.mute and cycle % 5 == 0:
            print ' Cycle', cycle, '  Current Score:', sf(pose)
        if args.pm:
            pm.apply(pose)
        if args.make_movie:
            pose.dump_pdb(dir_name + "/" + output_frame_base_filename +
                                                           str(cycle) + '.pdb')

    print ' Final Score:', sf(pose)
    sf.show(pose)  # TEMP
    if args.pm:
        pm.send_hbonds(pose)

    # Output results.
    pose.dump_pdb(output_pdb_filename)

    print 'Data written to:', output_pdb_filename


if __name__ == '__main__':
    # Parse arguments.
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('pdb_filename',
                        help='the filename of the PDB structure to be' +
                        ' evaluated')
    parser.add_argument('--mute', action='store_true',
                        help='flag to mute output for cycles')
    parser.add_argument('--make_movie', action='store_true',
                        help='flag to output a movie frame for each step')
    parser.add_argument('--pm', action='store_true',
                        help='flag to observe the protocol in PyMOL')
    parser.add_argument('--mm', action='store_true',
                        help='flag to use the molecular mechanics' +
                        ' score function')
    parser.add_argument('--lock_rings', action='store_true',
                        help='flag to prevent ring sampling from occurring')
    parser.add_argument('--n_cycles', type=int, default=100,
                        help='the number of Monte Carlo cycles')
    parser.add_argument('--kt', type=float, default=1.0,
                        help='the "temperature" to use for the MC cycles')
    parser.add_argument('--n_small_moves', type=int, default=5,
                        help='the number of small moves for each MC cycle')
    parser.add_argument('--n_shear_moves', type=int, default=5,
                        help='the number of shear moves for each MC cycle')
    parser.add_argument('--max_angle', type=float, default=10.0,
                        help='the maximum angle for a backbone move')
    parser.add_argument('--Hbond_weight', type=float, default=2.0,
                        help='the scoring weight for hydrogen bonds')
    parser.add_argument('--min_type', default='linmin',
                        choices=['linmin', 'dfpmin'],
                        help='the type of minimization')
    args = parser.parse_args()

    # Initialize Rosetta.
    init(extra_options='-include_sugars -read_pdb_link_records '
                       '-override_rsd_type_limit -mute all')

    # Create pose.
    print 'Loading pose from:', args.pdb_filename
    pose = pose_from_pdb(args.pdb_filename)

    # Set filenames.
    output_base_filename = args.pdb_filename[:-4] + '_refined'

    # Run protocol.
    refine_saccharide(pose, args, output_base_filename)