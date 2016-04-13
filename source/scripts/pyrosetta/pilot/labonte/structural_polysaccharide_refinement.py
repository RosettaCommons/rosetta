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
"""Brief:   This PyRosetta script refines an oligo- or polysaccharide chain.

Author:  Jason W. Labonte

"""
# Imports
from os import mkdir
from argparse import ArgumentParser

from rosetta import init, \
                    Pose, pose_from_pdb, \
                    get_fa_scorefxn, create_score_function, \
                    hbond_sr_bb, hbond_bb_sc, hbond_sc, hbond_lr_bb, \
                    fa_elec, fa_sol, fa_atr, fa_rep, fa_intra_rep, \
                    MoveMap, standard_packer_task, \
                    RotamerTrialsMover, SmallMover, ShearMover, MinMover, \
                    MonteCarlo, PyJobDistributor, \
                    PyMOL_Mover

from rosetta.core.chemical.carbohydrates import CarbohydrateInfo
from rosetta.protocols.simple_moves import RingConformationMover


# Global data
frame = 0  # for naming movie frame files


# Define methods.
def parse_arguments():
    """Use an ArgumentParser to parse the arguments."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('pdb_filename',
                        help='the filename of the PDB structure to be' +
                        ' evaluated')
    parser.add_argument('--output_folder', type=str,
                        default='output/refinement/',
                        help='the directory name where output files should ' +
                             'be saved')
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
    parser.add_argument('--n_decoys', type=int, default=1,
                        help='the number of docking decoys ' +
                             '(Setting to 0 simply pre-packs the structure.)')
    parser.add_argument('--kt', type=float, default=1.0,
                        help='the "temperature" to use for the MC cycles')
    parser.add_argument('--n_small_moves', type=int, default=5,
                        help='the number of small moves for each MC cycle')
    parser.add_argument('--n_shear_moves', type=int, default=5,
                        help='the number of shear moves for each MC cycle')
    parser.add_argument('--max_angle', type=float, default=10.0,
                        help='the maximum angle for a backbone move')
    parser.add_argument('--Hbond_mult', type=float, default=2.0,
                        help='the scoring weight multiplier for H-bonds')
    parser.add_argument('--elec_mult', type=float, default=1.0,
                        help='the scoring weight multiplier for electro' +
                             'static effects')
    parser.add_argument('--atr_mult', type=float, default=1.0,
                        help='the scoring weight multiplier for attractive ' +
                             'effects')
    parser.add_argument('--rep_mult', type=float, default=1.0,
                        help='the scoring weight multiplier for repulsive ' +
                             'effects')
    parser.add_argument('--sol_mult', type=float, default=1.0,
                        help='the scoring weight multiplier for solvent ' +
                             'effects')
    parser.add_argument('--min_type', default='linmin',
                        choices=['linmin', 'dfpmin'],
                        help='the type of minimization')
    return parser.parse_args()


def visualize(pose):
    if args.pm:
        pm.apply(pose)
    if args.make_movie:
        global frame
        path = args.output_folder + 'movie_frames/'
        if frame == 0:
            try:
                mkdir(path)
            except OSError:
                print 'Warning:',
                print 'Directory already exists;',
                print 'movie frame files will be overwritten.'
        frame += 1
        path += 'frame' + str(frame) + '.pdb'
        pose.dump_pdb(path)


def refine_saccharide(pose, args):
    """Perform packing, small, and shear moves."""
    print ' Initial Score:', sf(pose)

    # Set up BB Movers.
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(False)
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
        ring_min_mover = MinMover(mm_rings, sf, args.min_type, 0.001, True)

    # Set up Monte Carlo object.
    mc = MonteCarlo(pose, sf, args.kt)

    for cycle in range(1, args.n_cycles + 1):
        small_mover.apply(pose)
        shear_mover.apply(pose)
        visualize(pose)

        if not args.lock_rings:
            ring_flipper.apply(pose)
            ring_min_mover.apply(pose)
            visualize(pose)

        packer.apply(pose)
        min_mover.apply(pose)
        visualize(pose)

        mc.boltzmann(pose)
        if not args.mute and cycle % 5 == 0:
            print ' Cycle', cycle, '  Current Score:', sf(pose)
        visualize(pose)

    print ' Final Score:', sf(pose)
    sf.show(pose)  # TEMP
    #if args.pm:
    #    pm.send_hbonds(pose)


if __name__ == '__main__':
    # Parse arguments.
    args = parse_arguments()

    # Initialize Rosetta.
    print '\nInitializing Rosetta...'
    init(extra_options='-include_sugars '
                       '-write_pdb_link_records '
                       '-mute basic -mute numeric -mute utility '
                       '-mute core -mute protocols '
                       '-missing_density_to_jump '
                       '-default_max_cycles 200 '
                       '-default_repeats 30 '
                       '-set_weights cart_bonded 0.05 pro_close 0 '
                       '-minimize_bond_angles '
                       '-minimize_bond_lengths '
                       '-nonideal '
                       )

    # Create pose.
    print 'Loading pose from:', args.pdb_filename
    pose = pose_from_pdb(args.pdb_filename)

    # Set up ScoreFunction.
    if args.mm:
        sf = create_score_function('mm_std')
    else:
        sf = get_fa_scorefxn()
    sf.set_weight(hbond_sr_bb, sf.get_weight(hbond_sr_bb) * args.Hbond_mult)
    sf.set_weight(hbond_lr_bb, sf.get_weight(hbond_lr_bb) * args.Hbond_mult)
    sf.set_weight(hbond_bb_sc, sf.get_weight(hbond_bb_sc) * args.Hbond_mult)
    sf.set_weight(hbond_sc, sf.get_weight(hbond_sc) * args.Hbond_mult)
    sf.set_weight(fa_atr, sf.get_weight(fa_atr) * args.atr_mult)
    sf.set_weight(fa_intra_rep, sf.get_weight(fa_rep) * args.rep_mult)
    sf.set_weight(fa_sol, sf.get_weight(fa_sol) * args.sol_mult)
    if not args.mm:
        sf.set_weight(fa_elec, sf.get_weight(fa_elec) * args.elec_mult)

    # Set filenames.
    output_base_filename = args.pdb_filename[:-4] + '_refined'

    if args.pm:
        pose.pdb_info().name(output_base_filename)
        pm = PyMOL_Mover()

    visualize(pose)

    # Set up packer.
    pt = standard_packer_task(pose)
    pt.restrict_to_repacking()
    pt.or_include_current(True)

    packer = RotamerTrialsMover(sf, pt)

    # Prepack.
    packer.apply(pose)
    visualize(pose)
    print ' Score after prepacking:', sf(pose)
    sf.show(pose)  # TEMP

    # Fold the polysaccharide.
    print 'Folding/Refining structural polysaccharide with',
    print pose.total_residue(),
    #print pose.residue(2).carbohydrate_info().short_name(), 'residues...'
    print 'residues...'

    # Prepare job distributor.
    jd = PyJobDistributor(args.output_folder + output_base_filename,
                          args.n_decoys, sf)

    starting_pose = Pose()
    starting_pose.assign(pose)

    while not jd.job_complete:
        print ' Decoy', jd.current_num
        print '  Randomizing positions...'
        pose.assign(starting_pose)

        # Run protocol.
        print '  Refining...'
        refine_saccharide(pose, args)

        jd.output_decoy(pose)
        args.make_movie = False  # Make no more than one movie.