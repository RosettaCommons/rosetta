#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   beta_nonlocal.py
## @brief  Look for 'secondary structures' like sheets in beta AAs
## @author Andrew Watkins

if __name__ == '__main__':
    from pyrosetta import *
    from pyrosetta.rosetta import *

    init()

    # create score function
    score_fxn = core.scoring.get_score_function()
    score_fxn.set_weight(core.scoring.fa_rep, 0.3)
    score_fxn.set_weight(core.scoring.atom_pair_constraint, 1)
    score_fxn.set_weight(core.scoring.dihedral_constraint, 0.3)
    #score_fxn.set_weight(hbond_lr_bb, 50)
    #score_fxn.set_weight(hbond_sr_bb, 0.2)
    from rosetta.core.conformation import Residue, ResidueFactory

    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()

    pose = Pose()

    pose.append_residue_by_jump(Residue(rts.name_map("B3A:AcetylatedNtermProteinFull"), True), 1)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("NME"), True), True)

    pose.append_residue_by_jump(Residue(rts.name_map("B3A:AcetylatedNtermProteinFull"), True), 1, "", "", True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A"), True), True)
    pose.append_residue_by_bond(Residue(rts.name_map("B3A:MethylatedCtermProteinFull"), True), True)

    from rosetta.core.id import TorsionID, AtomID

    for ii in xrange(pose.size()):
        pose.set_torsion(TorsionID(ii + 1, core.id.BB, 4), 180)

    translate = protocols.rigid.RigidBodyTransMover(pose, 1)
    translate.step_size(4)
    translate.apply(pose)

    mm = core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_jump(True)
    for ii in xrange(pose.size()):
        mm.set(TorsionID(ii + 1, core.id.BB, 4), False)

    minM = protocols.simple_moves.MinMover(mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, True)
    minM.apply(pose)
    
    from rosetta.core.scoring.constraints import AtomPairConstraint, Constraint, MultiConstraint, AmbiguousConstraint
    from rosetta.core.scoring.func import HarmonicFunc
    
    ambig = AmbiguousConstraint()

    multi = MultiConstraint()
    for ii in xrange(1, 6, 2):
        jj = 12 - ii
        multi.add_individual_constraint(AtomPairConstraint(
            AtomID(pose.residue(ii).atom_index("O"), ii),
            AtomID(pose.residue(jj).atom_index("H"), jj),
            core.scoring.func.HarmonicFunc(1.9, 0.1)))
        multi.add_individual_constraint(AtomPairConstraint(
            AtomID(pose.residue(ii).atom_index("H"), ii),
            AtomID(pose.residue(jj).atom_index("O"), jj),
            core.scoring.func.HarmonicFunc(1.9, 0.1)))
    ambig.add_individual_constraint(multi)

    multi2 = MultiConstraint()
    for ii in xrange(1, 6, 2):
        jj = ii + 6
        multi2.add_individual_constraint(AtomPairConstraint(
            AtomID(pose.residue(ii).atom_index("O"), ii),
            AtomID(pose.residue(jj).atom_index("H"), jj),
            core.scoring.func.HarmonicFunc(1.9, 0.1)))
        if ii == 5: continue
        multi2.add_individual_constraint(AtomPairConstraint(
            AtomID(pose.residue(ii).atom_index("H"), ii),
            AtomID(pose.residue(jj+1).atom_index("O"), jj),
            core.scoring.func.HarmonicFunc(1.9, 0.1)))
    ambig.add_individual_constraint(multi)
    pose.add_constraint(ambig)

    pose.dump_pdb("out.pdb")
    for phi in xrange(-170, 190, 10):
        for tht in xrange(-170, 190, 10):
            for psi in xrange(-170, 190, 10):
                for ii in xrange(pose.size()):
                    if ii + 1 == 6: continue
                    from rosetta.core.id import TorsionID
                    pose.set_torsion(TorsionID(ii + 1, core.id.BB, 1), phi)
                    pose.set_torsion(TorsionID(ii + 1, core.id.BB, 2), tht)
                    pose.set_torsion(TorsionID(ii + 1, core.id.BB, 3), psi)

                minM.apply(pose)
                if pose.energies().total_energy() < -40:
                    print phi, tht, psi, pose.energies().total_energy()
                    pose.dump_scored_pdb("out_{}_{}_{}.pdb".format(phi, tht, psi), score_fxn)
