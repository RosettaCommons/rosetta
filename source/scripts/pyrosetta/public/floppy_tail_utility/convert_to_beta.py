#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

"""
Brief: Script to change the secondary structure of all n- or c-termini in a PDB.

Params: convert_to_beta.py –i <in.pdb> –s <residue #> –c –o <out.pdb>

Example: convert_to_beta.py –i 1hk9.pdb –s 65 –c –o 1hk9.extend.pdb --publication-specific

Remarks: Assumes that changes should be done to the n-/c-termini of all chains, but that may 
         not always be the case.

Author: Jeliazko Jeliazkov
"""
import pyrosetta
import rosetta
import argparse

def parse_args():
    # get start of beta structured terminus here
    # assume symmetry

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-i","--infile", type=str, help="input pdb file to adjust secondary structure", required=True)

    parser.add_argument('-s','--start-resn', type=int, help="residue from which to begin converting secondary structure, direction is specified by -c or -n.")

    parser.add_argument("-o","--outfile", type=str, help="output filename")

    parser.add_argument("--publication-specific", action='store_true') # use hard-coded values for angles to match xtal

    switch = parser.add_mutually_exclusive_group(required=True)
    switch.add_argument("-c","--c-termini",action='store_true')
    switch.add_argument("-n","--n-termini",action='store_true')

    args = parser.parse_args()
    return args

def n_chains(pose):
    """
    Return number / chains of a pose.
    """
    chains = set()
    for i in range(1,pose.size()+1):
        chains |= set(pose.pdb_info().chain(i))

    return chains

def chain_length(pose, chain):
    """
    Return length of a chain in pose.
    """
    chain_length = 0
    for i in range(1,pose.size()+1):
        if pose.pdb_info().chain(i) == chain:
            chain_length += 1

    return chain_length

if __name__ == '__main__':
    args = parse_args()
    print args

    # initialize Rosetta
    pyrosetta.init()

    # load pose
    pose = pyrosetta.pose_from_pdb(args.infile)

    # loop over chains and modify secondary structure
    chains = n_chains(pose)

    # updated jumps from n-terminal to somewhere else?
    # slide fold_tree jumps so n-terminal changes do not propegate
    if args.n_termini:
        ft = pose.fold_tree()

        for cn in range(len(chains)-1):
            ft.slide_jump(cn+1,args.start_resn+1,ft.jump_edge(cn+1).stop()+args.start_resn)
            pose.fold_tree(ft)

        print pose.fold_tree()

    for cn, chain in enumerate(sorted(chains)):

        print cn, chain, chain_length(pose,chain)

        if args.c_termini:

            start = pose.pdb_info().pdb2pose(chain,args.start_resn)
            stop = pose.pdb_info().pdb2pose(chain,chain_length(pose,chain))+1

            # don't change ft, default is good?

            for i in range(start,stop):

                pose.set_phi(i,-135)
                pose.set_psi(i,135)

                if args.publication_specific: # match 1HK9 @ 65-70
                    # hard-coded offset for 4NL2 that differs in numbering
                    #offset = 2 
                    offset = 0
                    if i == 64+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-60.88)
                        pose.set_psi(i,149.37)

                    if i == 65+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-82.0)
                        pose.set_psi(i,-10.2)

                    if i == 66+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-155.1)
                        pose.set_psi(i,160.5)

                    if i == 67+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-55.9)
                        pose.set_psi(i,142.1)

                    if i == 68+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-131.2)
                        pose.set_psi(i,107.5)

                    if i == 69+cn*chain_length(pose,chain)+offset:
                        pose.set_phi(i,-76.3)
                        pose.set_psi(i,160.0)
                    #if i == 70+cn*chain_length(pose,chain):
                    #    print pose.residue(i).name()
                    #    pose.set_phi(i,58.9)
                    #    pose.set_psi(i,87.2)

        elif args.n_termini: # start is higher than sotp

            start = pose.pdb_info().pdb2pose(chain,1)-1
            stop = pose.pdb_info().pdb2pose(chain,args.start_resn)

            for i in range(stop,start,-1):

                pose.set_phi(i,-135)
                pose.set_psi(i,135)

    if args.outfile:
        pose.dump_pdb(args.outfile)
    else:
        pose.dump_pdb(pose.pdb_info().name()+'.beta')
