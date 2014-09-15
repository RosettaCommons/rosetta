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

"""Brief:   This PyRosetta script performs a bound-bound dock between an
         oligosaccharide and a carbohydrate-binding protein.

Author:  Jason W. Labonte

"""

# Imports
from argparse import ArgumentParser
from os import remove, listdir
from math import exp
from random import random

from rosetta import init, Pose, pose_from_pdb, PyMOL_Mover, setup_foldtree, \
                    Vector1, get_fa_scorefxn, PyJobDistributor, MoveMap, \
                    create_score_function_ws_patch, DockMCMProtocol, \
                    FaDockingSlideIntoContact, MinMover, MonteCarlo, \
                    hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, fa_elec, \
                    fa_sol, fa_atr, fa_rep, standard_packer_task, \
                    TaskFactory, InitializeFromCommandline, IncludeCurrent, \
                    RestrictToRepacking, PackRotamersMover
from rosetta.protocols.rigid import RigidBodyPerturbMover, \
                                    RigidBodyRandomizeMover, \
                                    partner_upstream, partner_downstream
from rosetta.core.scoring import non_peptide_heavy_atom_RMSD


# Constants
JUMP_NUM = 1
STARTING_RAMP_DOWN_FACTOR = 3.25  # values used by Krishna
STARTING_RAMP_UP_FACTOR = 0.45455  # values used by Krishna


# Add custom MonteCarlo method.
def _distance_criterion(self, pose, last_distance, last_pose):
    current_distance = pose.jump(JUMP_NUM).get_translation().length
    distance_delta = current_distance - last_distance
    boltz_factor =  -distance_delta / self.temperature()
    probability = exp(min(40.0, max(-40.0, boltz_factor)))
    if probability < 1:
        #print '     Centers moved apart ( from', last_distance,
        #print 'to', current_distance, '):',
        if random() >= probability:
            #print 'rejecting'
            accepted = False
        else:
            #print 'accepting'
            accepted = True
    else:
        #print '     Centers moved closer ( from', last_distance,
        #print 'to', current_distance, '): accepting'
        accepted = True

    if not accepted:
        pose.assign(last_pose)
        self.set_last_accepted_pose(last_pose)
    else:
        self.set_last_accepted_pose(pose)
    return accepted

MonteCarlo.distance_criterion = _distance_criterion


# Define methods.
def parse_arguments():
    """Use an ArgumentParser to parse the arguments."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('pdb_filename1',
                        help='the filename of one of the PDB structures to ' +
                             'be docked')
    parser.add_argument('pdb_filename2',
                        help='the filename of one of the PDB structures to ' +
                             'be docked')
    parser.add_argument('--output_folder', type=str,
                        default='output/docking/global/',
                        help='the directory name where output files should ' +
                             'be saved')
    parser.add_argument('--local', action='store_true',
                        help='perform local docking')
    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite existing output files')
    parser.add_argument('--prepack', action='store_true',
                        help='prepack the input structures')
    parser.add_argument('--ref',
                        help='the filename of the PDB structure to be used ' +
                             'as a reference')
    parser.add_argument('--mute', action='store_true',
                        help='flag to mute output for cycles')
    parser.add_argument('--make_movie', action='store_true',
                        help='flag to output a movie frame for each step')
    parser.add_argument('--pm', action='store_true',
                        help='flag to observe the protocol in PyMOL')
    parser.add_argument('--mm', action='store_true',
                        help='flag to use the molecular mechanics' +
                        ' score function')
    parser.add_argument('--n_cycles', type=int, default=100,
                        help='the number of Monte Carlo refinement cycles')
    parser.add_argument('--kt', type=float, default=0.8,
                    help='the "temperature" to use for the MC cycles')
    parser.add_argument('--n_decoys', type=int, default=10,
                        help='the number of docking decoys ' +
                        '(Setting to 0 simply pre-packs the structure.)')
    parser.add_argument('--rot', type=float, default=2.0,
                        help='the mean rotation (in degrees) for rigid-body ' +
                             'moves')
    parser.add_argument('--trans', type=float, default=0.5,
                        help='the mean translation (in Angstroms) for rigid-' +
                             'body moves')
    parser.add_argument('--Hbond_mult', type=float, default=1.0,
                        help='the scoring weight multiplier for hydrogen ' +
                             'bonds')
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
    return parser.parse_args()


def clrdir(directory):
    """Delete all files in the given directory."""
    for f in listdir(directory):
        remove(directory + f)


def merge_pdb_files(filename1, filename2, output_folder):
    """Merges two PDB files and returns the path of the new file."""
    if not filename1.endswith('.pdb') or not filename2.endswith('.pdb'):
        exit('Input files must have the ".pdb" file extension.')
    combined_filename = output_folder + filename1[:-4] + "-" + filename2
    with open(combined_filename, "w") as new_file:
        with open(filename1) as file1:
            new_file.write(file1.read())  # copy file1
        with open(filename2) as file2:
            new_file.write(file2.read())  # append file2
    return combined_filename


def visualize(pose):
    if args.pm:
        pm.apply(pose)
    if args.make_movie:
        pass  # need a global frame count


def ramp_score_weight(score_function, method, target, fraction_completion):
    #print 'fraction completion:', fraction_completion, ":",  # DEBUG
    current_weight = score_function.get_weight(method)
    #print 'current:', current_weight,
    #print 'target:', target,
    if current_weight < target:
        factor_range_size = 1 - STARTING_RAMP_UP_FACTOR
        factor = STARTING_RAMP_UP_FACTOR + \
                                        fraction_completion * factor_range_size
    elif current_weight > target:
        factor_range_size = STARTING_RAMP_DOWN_FACTOR - 1
        factor = STARTING_RAMP_DOWN_FACTOR - \
                                        fraction_completion * factor_range_size
    else:
        factor = 1
    new_weight = target * factor
    #print 'new:', new_weight
    score_function.set_weight(method, new_weight)


if __name__ == '__main__':
    # Parse arguments.
    args = parse_arguments()

    if args.overwrite:
        print '\nClearing out directory:', args.output_folder + '...',
        clrdir(args.output_folder)
        print 'Directory cleared.'

    # Merge input files.
    print '\nGenerating starting structure file...',
    new_filename = merge_pdb_files(args.pdb_filename1,
                                   args.pdb_filename2,
                                   args.output_folder)
    print ' Success!  Created new file:', new_filename

    # Initialize Rosetta.
    print '\nInitializing Rosetta...'
    init(extra_options='-include_sugars -include_lipids ' +
                       '-read_pdb_link_records -write_pdb_link_records ' +
                       '-override_rsd_type_limit -mute basic -mute numeric ' +
                       '-mute utility -mute core -mute protocols')

    # Create pose.
    print '\nGenerating starting pose...'
    starting_pose = pose_from_pdb(new_filename)

    # Display stating pose.
    if args.pm:
        starting_pose.pdb_info().name(args.pdb_filename1[:-4] + "-" + \
                                                            args.pdb_filename2)
        pm = PyMOL_Mover()
    visualize(starting_pose)

    # Prepare foldtree
    chain_endings = starting_pose.conformation().chain_endings()
    partners = starting_pose.pdb_info().chain(chain_endings[JUMP_NUM]) + \
              "_" + starting_pose.pdb_info().chain(chain_endings[JUMP_NUM] + 1)

    setup_foldtree(starting_pose, partners, Vector1([JUMP_NUM]))

    # Print some information about the starting pose.
    print "", starting_pose.fold_tree(),
    print 'Ligand center is',
    print starting_pose.jump(JUMP_NUM).get_translation().length,
    print 'Angstroms from protein center.'

    # Prepare scoring function.
    #sf = get_fa_scorefxn()
    sf = create_score_function_ws_patch("talaris2013", "docking")
    sf.set_weight(hbond_sr_bb, sf.get_weight(hbond_sr_bb) * args.Hbond_mult)
    sf.set_weight(hbond_lr_bb, sf.get_weight(hbond_lr_bb) * args.Hbond_mult)
    sf.set_weight(hbond_bb_sc, sf.get_weight(hbond_bb_sc) * args.Hbond_mult)
    sf.set_weight(hbond_sc, sf.get_weight(hbond_sc) * args.Hbond_mult)
    sf.set_weight(fa_elec, sf.get_weight(fa_elec) * args.elec_mult)
    sf.set_weight(fa_sol, sf.get_weight(fa_sol) * args.sol_mult)
    target_atr = sf.get_weight(fa_atr) * args.atr_mult
    target_rep = sf.get_weight(fa_rep) * args.rep_mult
    sf.set_weight(fa_atr, target_atr)
    sf.set_weight(fa_rep, target_rep)
    print ' Starting Score:'
    sf.show(starting_pose)

    if args.prepack:
        # Prepare PackerTask
        pt = standard_packer_task(starting_pose)
        pt.restrict_to_repacking()
        pt.or_include_current(True)

        # Prepare the packer.
        pre_packer = PackRotamersMover(sf, pt)

        # Pre-pack the pose.
        print '\nPre-packing the pose (slow)...'
        pre_packer.apply(starting_pose)
        print ' Success!'
        print ' Pre-packed score:',
        sf.show(starting_pose)
        starting_pose.dump_pdb(new_filename[:-4] + 'pre-packed.pdb')
        visualize(starting_pose)

    if not args.n_cycles:
        exit()

    # Prepare other Movers.
    slider = FaDockingSlideIntoContact(JUMP_NUM)
    perturber = RigidBodyPerturbMover(JUMP_NUM, args.rot, args.trans)
    randomizerA = RigidBodyRandomizeMover(starting_pose, JUMP_NUM,
                                          partner_upstream)
    randomizerB = RigidBodyRandomizeMover(starting_pose, JUMP_NUM,
                                          partner_downstream)
    mm = MoveMap()
    mm.set_jump(JUMP_NUM, True)
    jump_minimizer = MinMover()
    jump_minimizer.movemap(mm)
    jump_minimizer.score_function(sf)
    #dock_hires = DockMCMProtocol()
    #dock_hires.set_scorefxn(sf)
    #dock_hires.set_partners(partners)

    # Prepare job distributor.
    jd = PyJobDistributor(new_filename[:-4], args.n_decoys, sf)

    if args.ref:
        if not args.ref.endswith('.pdb'):
            exit('Reference file must have the ".pdb" file extension.')
        ref_pose = pose_from_pdb(args.ref)
        #jd.native_pose = ref_pose

    # Begin docking protocol.
    print "\nDocking..."
    pose = Pose()  # working pose
    last_pose = Pose()  # needed because I have not modified C++ code
    while not jd.job_complete:
        print ' Decoy', jd.current_num
        print '  Randomizing positions...'
        pose.assign(starting_pose)
        print 'JUMP LENGTH:', pose.jump(JUMP_NUM).get_translation().length  # DEBUG
        if not args.local:
            randomizerA.apply(pose)
        randomizerB.apply(pose)
        visualize(pose)
        slider.apply(pose)
        visualize(pose)
        print '   Ligand center is',
        print pose.jump(JUMP_NUM).get_translation().length,
        print 'Angstroms from protein center.'

        print '  Refining...'
        # Set beginning values for weights to ramp.
        sf.set_weight(fa_atr, target_atr * STARTING_RAMP_DOWN_FACTOR)
        sf.set_weight(fa_rep, target_atr * STARTING_RAMP_UP_FACTOR)
        mc = MonteCarlo(pose, sf, args.kt)
        for cycle in range(1, args.n_cycles + 1):
            if cycle % (args.n_cycles/10) == 0:  # Ramp every ~10% of n_cycles
                fraction = float(cycle) / args.n_cycles
                ramp_score_weight(sf, fa_atr, target_atr, fraction)
                ramp_score_weight(sf, fa_rep, target_rep, fraction)
                mc.reset(pose)
            # Save information needed for distance criterion
            #old_score = sf(pose)
            distance = pose.jump(JUMP_NUM).get_translation().length
            last_pose.assign(pose)
            perturber.apply(pose)
            visualize(pose)
            #jump_minimizer.apply(pose)  # Seg-faults; why?
            #visualize(pose)
            #print '    Score changed from', old_score, 'to', sf(pose), ':',
            if mc.boltzmann(pose):
                # Only perform a distance acceptance check if the new pose is
                # accepted.
                #print 'accepting pending distance check:'
                mc.distance_criterion(pose, distance, last_pose)
            #else:
                #print 'rejecting'
            if not args.mute and cycle % 5 == 0:
                print '   Cycle', cycle, '  Current Score:', sf(pose)
            visualize(pose)

        print '  Final Score for decoy', jd.current_num, ':', sf(pose)
        #sf.show(pose)  # TEMP
        if args.ref:
            rmsd = non_peptide_heavy_atom_RMSD(pose, ref_pose)
            jd.additional_decoy_info = 'ligand_rmsd: ' + str(round(rmsd, 2))
        jd.output_decoy(pose)
