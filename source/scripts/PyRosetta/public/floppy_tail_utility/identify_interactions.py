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
Brief: Computes 2-body energies between two sets of residues for PDBs in a dir.

Params: identify_interactions.py --set1 <res#-res#> --set2 <res#,...,res#> -n_procs <# of processes to run> </path/to/pdbs/>

Example: identify_interactions.py --pdb_numbering --n_procs 14 --set1 66-82 --set2 1-65 --manual_suffix .pdb.gz decoys

Remarks: One can specify either Rosetta numbering (default) or a funky pdb numbering,
         where the sets are expended to include all chains with those residue numbers.

Author: Jeliazko Jeliazkov
"""

import pyrosetta
import rosetta
import argparse
import os
import sys
import itertools
from multiprocessing import Process

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("pdb_path", type=str, help="directory full of pdbs to analyze")

    parser.add_argument("--pdb_numbering", action='store_true',
                        help='Output chain + pdb number rather than pose number')
    parser.add_argument("--manual_suffix", type=str, default=".pdb",
                        help='manually specify suffix if not *.pdb')
    parser.add_argument("--output_file", type=str, default="2b_energies")
    parser.add_argument("--n_procs", default=1, type=int)

    parser.add_argument('--set1', type=str, help='PDB select (range 1-6 or csv 1,2,5)', required=True)
    parser.add_argument('--set2', type=str, help='PDB select (range 7-21 or csv 7,15)', required=True)

    # idea here is if the thresh is -1.0, we take all interactions with E < -1.0
    parser.add_argument('--threshold', type=float, help='Set cutoff for defining an intearction. Absolute values will be compared.', default=-1.0)

    args = parser.parse_args()

    return args

def parse_commandline_residues(cmd_residues):

    if cmd_residues.find('-') != -1:
        residues = [int(x) for x in cmd_residues.split('-')]
        residues = range(residues[0],residues[1]+1)
    elif cmd_residues.find(',') != -1:
        residues = [int(x) for x in cmd_residues.split(',')]
    else:
        try:
            residues = int(cmd_residues)
        except ValueError:
            print 'Could not parse command line residues!'
            sys.exit(1)

    return residues

def renumber_set(pdbfile,residue_set):
    """
    Take in residue numbers and expand to all chains.
    Mostly for Hfq stuff. Will need further refinining.
    16,17,19 should become 16A, 17A, 19A, 16B, 17B, 19B ...
    in pose numbering for all chains
    """
    renumbered_set = []

    pose = pyrosetta.pose_from_pdb(pdbfile)
    chains = get_chains(pose)

    for chain in chains:
        for res in residue_set:
            pose_num = pose.pdb_info().pdb2pose(chain,res)
            if pose_num != 0: renumbered_set.append(pose_num)

    return renumbered_set

def get_chains(pose):
    chains = set()
    for i in range(1,pose.size()+1):
        chain = pose.pdb_info().chain(i)
        if chain not in chains: chains = chains | set(chain)
    return chains


def chunkify(l,n):
    """
    Split list l into n chunks.
    """
    # minimum chunk length (last chunk may be longer)
    cl = len(l)/n
    chunks = []
    for i in range(n-1):
        chunks.append(l[i*cl:(i+1)*cl])
    chunks.append(l[(n-1)*cl:])
    return chunks

def calc_interaction_energies(pose, sfxn, scores_types, res1, res2):
    """
    Using Rosetta's energy function, this function calculates residue pairwise
    interactions in a protein and then reports the types & scores.
    """
    # init emap for storing values
    emap = rosetta.core.scoring.EMapVector()

    # vector of scores, order is the same?
    resultant_scores = []

    # might need to run eval_long_range_twobody_energies to get dslf?
    sfxn.eval_ci_2b(pose.residue(res1), pose.residue(res2), pose, emap)
    sfxn.eval_cd_2b(pose.residue(res1), pose.residue(res2), pose, emap)

    for sf_type in scores_types:
        resultant_scores.append(emap[sf_type]) # extract value of score type

    # personally, I think this is silly, but it needed; see ResidueScoreFeatures
    emap.clear()

    return resultant_scores

def calc_all_interaction_energies(pdbs, resnum_range1, resnum_range2, threshold, out_f, proc_id):
    # assume default score function, but this should be editable in the future
    sfxn = pyrosetta.get_fa_scorefxn()

    # double looping in list comprehensions is insane... but very efficient since python optimizes line by line
    # see below for more info http://stackoverflow.com/questions/3471999/how-do-i-merge-two-lists-into-a-single-list
    score_types = [item for pair in zip(sfxn.ci_2b_types(),sfxn.cd_2b_types()) for item in pair if sfxn.has_nonzero_weight(item)]
    score_types.append(rosetta.core.scoring.fa_elec) # manual hack
    score_types_str = [str(x) for x in score_types]

    # fix this to be editable later
    with open('{0}.{1}.out'.format(out_f,proc_id),'w') as mf:

        # write total number of PDBs to be analyzed
        mf.write('# unique PDBs in: {0}\n'.format(len(pdbs)))

        # write legend
        mf.write('pdb_name,res1,res2')
        for score_type in score_types_str:
            mf.write(',{0}'.format(score_type))
        mf.write(',total_E\n')

        # loop over all pdbs and get pairwise energies
        for pdb in pdbs:

            pose = pyrosetta.pose_from_pdb(pdb)
            total_score = sfxn(pose) # update energies
            pdb_name = pose.pdb_info().name()

            for res1 in resnum_range1:
                for res2 in resnum_range2:
                    # ignore self interactions and immediate neighbors
                    if res1 == res2 or res1+1 == res2 or res1-1 == res2 \
                    or res1+2 == res2 or res1-2 == res2:
                        continue
                    #if res1 in resnum_range2 or res2 in resnum_range1:
                    #    continue

                    scores = calc_interaction_energies(pose,sfxn,score_types,res1,res2)

                    # check if any single energy beats the threshold
                    for score in scores:
                        if score < threshold: # good!
                            break
                    else: # loop terminated due to list exhaustion
                        continue # no score less than threshold, no interactions

                    # there is an interaction!
                    mf.write('{0},{1},{2}'.format(pdb_name, res1, res2))

                    for score in scores:
                        mf.write(',{0}'.format(score))

                    mf.write(',{0}\n'.format(total_score))

    return

if __name__=='__main__':
    uargs = parse_args()

    # Expand residue selections into list
    set1 = parse_commandline_residues(uargs.set1)
    set2 = parse_commandline_residues(uargs.set2)

    # init Rosetta
    pyrosetta.init()

    pdbs = [os.path.join(uargs.pdb_path,x) for x in os.listdir(uargs.pdb_path) if x.endswith(uargs.manual_suffix)]

    if uargs.pdb_numbering:  # will expand selection for each chain!
        set1 = renumber_set(pdbs[0],set1)
        set2 = renumber_set(pdbs[0],set2)

    print 'Using residue set {0} to identify interactions with set {1}'.format(set1,set2)

    # calc_all_interaction_energies goes here
    # split pdb list into n chunks
    #pdb_chunks = chunkify(args.pdb_in,args.n_procs)
    pdb_chunks = chunkify(pdbs,uargs.n_procs)

    for i in range(uargs.n_procs):
        p = Process(target=calc_all_interaction_energies, args=(pdb_chunks[i], set1, set2, uargs.threshold, uargs.output_file, i))
        p.start()

    p.join()

    # merge files produced by output... assume all files contain identical
    # header and only keep the header from the first file.
    output_files = [x for x in os.listdir("./") if x.startswith(uargs.output_file)]
    print "Combining files ", output_files
    header = False
    concat_file = open(uargs.output_file+".out","w")

    for ofile in output_files:
        with open(ofile,"r") as mf:
            # first time, copy over header
            if not header:
                mf.readline()
                concat_file.write(mf.readline())
                header = True
            else: # skip header ( 2 lines )
                mf.readline()
                mf.readline()

            # always loop over lines copying
            for line in mf:
                concat_file.write(line)
        os.remove(ofile)
