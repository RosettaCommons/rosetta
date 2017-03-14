#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Author: Jeliazko Jeliazkov

"""
Brief: Script which appends or prepends residues in an extended conformation.

Params: extend_terminus.py –c <chain> –o <output.pdb> –p/-a <input.pdb> <AA sequence>

Example: extend_terminus.py –c A –o 1hk9.Ap.pdb –p 1hk9.clean.pdb MAKGQ

Author: Jeliazko Jeliazkov
"""
import pyrosetta
import rosetta
import argparse
import sys
import os

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("filename", type=str, help="input pdb file to prepend/append")
    parser.add_argument("sequence", type=str, help="sequence to prepend or append")

    parser.add_argument("-c","--chain", type=str, help="chain to prepend or append")

    parser.add_argument("-o","--outfile", type=str, help="output filename")

    switch = parser.add_mutually_exclusive_group(required=True)
    switch.add_argument("-p","--prepend",action='store_true')
    switch.add_argument("-a","--append",action='store_true')

    args = parser.parse_args()

    return args

def find_chain_termini(pose,chain):
    """
    Helper function to identify the pose numbering of the chain start/end.
    Not very efficiently coded.
    """

    start = pose.size()+1
    end = 0

    for i in range(1,pose.size()+1):
        if pose.pdb_info().chain(i) == chain:
            if i > end:
                end = i
            if i < start:
                start = i

    # requires converting chain from char to int representation
    # alternatively use pose.conformation().chain_begin()
    # alternatively use pose.conformation().chain_end()

    return start, end

if __name__ == "__main__":

    args = parse_args()

    pyrosetta.init()

    # take user sequence, but flank with dummy residues so restypes are correct
    seq_in = args.sequence

    # generate pose from sequence
    #dummy_pose = rosetta.pose_from_sequence(seq_in)

    # load in working pose
    pose = pyrosetta.pose_from_pdb(args.filename)

    # test chain finding function
    if args.chain and len(args.chain)==1:
        chain_nterm, chain_cterm = find_chain_termini(pose,args.chain)
        # first arg is the start, second arg is the end
    else:
        sys.exit("No chain specified!")

    # let's test some things
    print "Loaded pose with {0} residues!".format(pose.size())

    # get type set to generate residues
    res_type_set = pose.residue_type_set_for_pose()

    # replace current terminal residue, with non-terminal variant
    # then extend terminus
    if args.append:
        for i in range(len(seq_in)):
            res_type = res_type_set.get_representative_type_name1(seq_in[i])
            res = rosetta.core.conformation.ResidueFactory.create_residue(res_type)
            pose.conformation().safely_append_polymer_residue_after_seqpos(res, chain_cterm + i, True)

    elif args.prepend:
        for i in range(len(seq_in)):
            res_type = res_type_set.get_representative_type_name1(seq_in[::-1][i])
            res = rosetta.core.conformation.ResidueFactory.create_residue(res_type)
            pose.conformation().safely_prepend_polymer_residue_before_seqpos(res, chain_nterm, True)

    # let's test some things
    print "Extended pose has {0} residues!".format(pose.size())

    # iterate over new residues and set bb angles, but don't move "old" residues
    movemap = pyrosetta.MoveMap()

    # regions to repack / extend are the tail
    if args.append: 
        movemap.set_bb_true_range(chain_cterm+1, chain_cterm+len(seq_in))
        movemap.set_chi_true_range(chain_cterm+1, chain_cterm+len(seq_in))

    elif args.prepend:
        # reverse fold tree to not perturb input
        ft = pose.fold_tree()
        ft.reorder(chain_nterm+len(seq_in))
        pose.fold_tree(ft)
        movemap.set_bb_true_range(chain_nterm,chain_nterm+len(seq_in)-1)
        movemap.set_chi_true_range(chain_nterm,chain_nterm+len(seq_in)-1)

    print movemap

    # let's test some things
    print "Setting ideal beta strand backbone angles for..."

    for i in range(1,pose.size()+1,1):
        if movemap.get_bb(i):
            print i, pose.pdb_info().pose2pdb(i)
            pose.set_phi(i,-135)
            pose.set_psi(i, 135)
            pose.set_omega(i,180)

    # score and repack
    sfxn = pyrosetta.get_fa_scorefxn()

    task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    task.restrict_to_repacking()

    pack_rotamers = rosetta.protocols.simple_moves.PackRotamersMover(sfxn, task)
    pack_rotamers.apply(pose)

    if args.outfile:
        outf = os.path.join(os.getcwd(),os.path.abspath(args.outfile))
    else:
        outf = os.path.join(os.getcwd(),"extended.pdb")

    pose.dump_pdb(outf)
