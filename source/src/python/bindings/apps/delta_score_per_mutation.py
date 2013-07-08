# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
"""Brief:   This PyRosetta script mutates every residue to each of the 20
      canonical amino acids and calculates the change in Rosetta score.

Authors:  Jason W. Labonte & Michael Pacella

"""
# Imports
import sys
import argparse
from rosetta import *


# Constants
PACK_RADIUS = 10.0
AAs = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")


def mutate_residue(pose, mutant_position, mutant_aa, pack_radius=0.0,
                                                             pack_scorefxn=''):
    """
    Replace the residue at <mutant_position> in <pose> with <mutant_aa> and re-
    pack any residues within <pack_radius> angstroms of the mutating residue's
    center (nbr_atom) using <pack_scorefxn>
    Note: <mutant_aa> is the single letter name for the desired ResidueType

    This function was written by Evan Baugh and modified.

    """
    if not pose.is_fullatom():
        IOError('mutate_residue() only works with full-atom poses.')

    test_pose = Pose()
    test_pose.assign(pose)

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = get_fa_scorefxn()

    task = standard_packer_task(test_pose)
    task.or_include_current(True)

    # A vector1 of booleans (a specific object) is needed for specifying the
    # mutation.  This demonstrates another more direct method of setting
    # PackerTask options for design.
    aa_bool = rosetta.utility.vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes).
    # The numbers 1-20 correspond individually to the 20 proteogenic amino
    # acids.  aa_from_oneletter_code() returns the integer representation of an
    # amino acid from its one letter code

    # Convert mutant_aa to its integer representation.
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # The mutation is performed by using a PackerTask with only the mutant
    # amino acid available during design.  To do this, we construct a vector1
    # of booleans indicating which amino acid (by its numerical designation;
    # see above) to allow.
    for i in range(1, 20 + 1):
        # In Python, logical expression are evaluated with priority, thus the
        # line below appends to aa_bool the truth (True or False) of the
        # statement i == mutant_aa.
        aa_bool.append(i == mutant_aa)

    # Modify the mutating residue's assignment in the PackerTask using the
    # vector1 of booleans across the proteogenic amino acids.
    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(
                                                                       aa_bool)

    # Prevent residues from packing by setting the per-residue "options" of
    # the PackerTask.
    center = pose.residue(mutant_position).nbr_atom_xyz()
    for i in range(1, pose.total_residue() + 1):
        # Only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                         test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:
            task.nonconst_residue_task(i).prevent_repacking()

    # Apply the mutation and pack nearby residues.
    packer = PackRotamersMover(pack_scorefxn, task)
    packer.apply(test_pose)

    return test_pose


# Parse arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('pdb_filename',
                    help='the filename of the PDB structure to be evaluated')
parser.add_argument('-m', '--minimize', action='store_true',
                    help='flag to perform minimization after each mutation')
args = parser.parse_args()

if __name__ == '__main__':
    rosetta.init(extra_options='-mute basic -mute core')

    data = ['Variant,Rosetta Score,"delta-delta-G"\n']

    # Load pdb file.
    initial_pose = pose_from_pdb(args.pdb_filename)

    # Set up ScoreFunction
    sf = get_fa_scorefxn()

    # Set up MoveMap
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    # Pack and minimize initial pose to remove clashes.
    pre_pre_packing_score = sf(initial_pose)

    task = standard_packer_task(initial_pose)
    task.restrict_to_repacking()
    task.or_include_current(True)
    pack_rotamers_mover = RotamerTrialsMover(sf, task)
    pack_rotamers_mover.apply(initial_pose)

    min_mover = MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(sf)
    min_mover.min_type('linmin')
    if args.minimize:
        min_mover.apply(initial_pose)

    post_pre_packing_score = sf(initial_pose)

    print
    print 'Reference Protein:', args.pdb_filename
    print '  Score:'
    print '    Before pre-packing:', pre_pre_packing_score
    print '    After pre-packing:', post_pre_packing_score
    print

    data.append('WT,' + str(post_pre_packing_score) + ',0.0\n')

    # Loop over all residues in the protein.
    print 'Making mutations and scoring poses. ',
    print '(This will take a while; please wait....)'

    n_res = initial_pose.total_residue()
    for seq_pos in range(1, n_res + 1):
        res = initial_pose.residue(seq_pos)
        for AA in AAs:
            if res.name1() != AA:
                variant_name = \
                        res.name1() + \
                        str(initial_pose.pdb_info().number(seq_pos)) + \
                        AA

                # Check for disulfide special case.
                if res.name() == 'CYD':
                    disulfide_partner = res.residue_connection_partner(
                                                   res.n_residue_connections())
                    temp_pose = Pose()
                    temp_pose.assign(initial_pose)
                    # (Packing causes seg fault if current CYS residue is not
                    # also converted before mutating.)
                    change_cys_state(seq_pos, 'CYS',
                                     temp_pose.conformation())
                    change_cys_state(disulfide_partner, 'CYS',
                                     temp_pose.conformation())
                    # Mutate protein.
                    mutant_pose = mutate_residue(temp_pose, seq_pos, AA,
                                                 PACK_RADIUS, sf)
                else:
                    # Mutate protein.
                    mutant_pose = mutate_residue(initial_pose, seq_pos, AA,
                                                 PACK_RADIUS, sf)

                # Minimize.
                if args.minimize:
                    min_mover.apply(mutant_pose)

                # Score.
                variant_score = sf(mutant_pose)

                data.append(variant_name + "," + \
                            str(variant_score) + "," + \
                            str(variant_score - post_pre_packing_score) + "\n")
        print '  All variants for 1st', seq_pos,
        print 'residues (of', n_res,
        print 'residues) made and scored.\r',
        sys.stdout.flush()
    print '\nMutations and scoring complete.'

    # Output results.
    data_filename = args.pdb_filename[:-4] + '_variant_scores.csv'
    with open(data_filename, "w") as f:
        f.writelines(data)

    print 'Data written to:', data_filename