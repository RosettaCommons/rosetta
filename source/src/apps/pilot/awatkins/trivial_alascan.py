#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   trivial_alascan.py
## @brief  An alascan that does no sampling
## @author Andrew Watkins

def mutate_residue_to_ala(pose, ii, score_fxn):

    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()
    res = pose.residue(ii)

    variant_types = utility.vector1_std_string(res.type().properties().get_list_of_variants())
    ala_type = rosetta.core.chemical.ResidueTypeFinder(rts).name3("ALA").variants(variant_types).get_representative_type()

    new_res = rosetta.core.conformation.ResidueFactory.create_residue(ala_type, res, pose.conformation())
    core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(res, new_res, pose.conformation())
    pose.conformation().replace_residue(ii, new_res, False)
    return score_fxn(pose)

def alascan(pose):
    score_fxn = core.scoring.get_score_function()
    core.scoring.constraints.add_fa_constraints_from_cmdline_to_scorefxn(score_fxn)
    score_fxn.set_weight(core.scoring.ref, 0)
    score_fxn.set_weight(core.scoring.unfolded, 0)

    wt_score = (score_fxn)(pose)
    print "Wildtype scores", wt_score

    # Don't even bother to calculate the interface yet
    for ii in xrange(pose.size()):

        pose_copy = Pose(pose)
        mut_score = mutate_residue_to_ala(pose_copy, ii + 1, score_fxn)

        if mut_score - wt_score > 2:
            print pose.pdb_info().pose2pdb(ii + 1), (mut_score - wt_score)
            


if __name__ == '__main__':
    from pyrosetta import *
    from pyrosetta.rosetta import *
    init()

    import sys
    pose = core.import_pose.pose_from_file(sys.argv[1])
    # pose on command line
    alascan(pose)

