#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=True:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   orn_macrocycle.py
## @brief  Create a macrocycle using ornithine residues
## @author Andrew Watkins


def add_orn_cst(pose, res1, res2):
    from rosetta.core.scoring.constraints import DihedralConstraint, AngleConstraint, AtomPairConstraint
    from rosetta.core.scoring.func import CircularHarmonicFunc, HarmonicFunc

    pose.add_constraint(
        DihedralConstraint( 
            AtomID(pose.residue(res1).atom_index("CA"), res1),
            AtomID(pose.residue(res1).atom_index("C" ), res1),
            AtomID(pose.residue(res2).atom_index("NE"), res2),
            AtomID(pose.residue(res2).atom_index("CD"), res2),
            CircularHarmonicFunc(3.14159, 0.02)))
    
    pose.add_constraint(
        DihedralConstraint( 
            AtomID(pose.residue(res1).atom_index("O"), res1),
            AtomID(pose.residue(res1).atom_index("C" ), res1),
            AtomID(pose.residue(res2).atom_index("NE"), res2),
            AtomID(pose.residue(res2).atom_index("1HE"), res2),
            CircularHarmonicFunc(3.14159, 0.02)))
    
    pose.add_constraint(
        DihedralConstraint( 
            AtomID(pose.residue(res1).atom_index("O"), res1),
            AtomID(pose.residue(res1).atom_index("C" ), res1),
            AtomID(pose.residue(res2).atom_index("CA"), res2),
            AtomID(pose.residue(res2).atom_index("NE"), res2),
            CircularHarmonicFunc(3.14159, 0.02)))

    pose.add_constraint(
        AngleConstraint( 
            AtomID(pose.residue(res1).atom_index("O"), res1),
            AtomID(pose.residue(res1).atom_index("C" ), res1),
            AtomID(pose.residue(res2).atom_index("NE"), res2),
            CircularHarmonicFunc(3.14159*2.0/3.0, 0.02)))


    pose.add_constraint(
        AtomPairConstraint( 
            AtomID(pose.residue(res1).atom_index("C" ), res1),
            AtomID(pose.residue(res2).atom_index("NE"), res2),
            HarmonicFunc(1.33, 0.02)))


if __name__ == '__main__':
    from pyrosetta import *
    from pyrosetta.rosetta import *

    from rosetta.core.conformation import Residue, ResidueFactory

    init()
    pose = Pose()
    score_fxn = core.scoring.get_score_function()
    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()
    ala_type = rts.name_map("ALA")
    orn_type = rts.name_map("ORN:NtermProteinFull:N-conjugated")
    ala = Residue(ala_type, True)
    orn = Residue(orn_type, True)

    pose.append_residue_by_jump(orn, 1)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_atoms(orn, True, "NE", 6, "C") 
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)
    pose.append_residue_by_bond(ala, True)

    # Hardcoded connection IDs; never do this! This just illustrates the use case.
    # Note that the NE connections are still #2 because they're ntermproteinfull.

    new_orn = ResidueFactory.create_residue(orn_type, pose.residue(1), pose.conformation())
    core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms(pose.residue(1), new_orn, pose.conformation())
    new_orn.residue_connection_partner( 2, 12, 2 )
    pose.conformation().replace_residue( 1, new_orn, False )
    new_ala = ResidueFactory.create_residue( ala_type, pose.residue( 12 ), pose.conformation() )
    core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( 12 ), new_ala, pose.conformation() )
    new_ala.residue_connection_partner( 2,  1, 2 )
    pose.conformation().replace_residue( 12, new_ala, False )
    
    pose.conformation().declare_chemical_bond( 12, "C", 1, "NE" )

    new_orn2 = ResidueFactory.create_residue( orn_type, pose.residue( 7 ), pose.conformation() )
    core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( 7 ), new_orn2, pose.conformation() )
    new_orn2.residue_connection_partner( 2,  6, 2 )
    pose.conformation().replace_residue( 7, new_orn2, False )
    
    new_ala2 = ResidueFactory.create_residue( ala_type, pose.residue(  6 ), pose.conformation() )
    core.conformation.copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue(  6 ), new_ala2, pose.conformation() )
    new_ala2.residue_connection_partner( 2,  7, 2 )
    pose.conformation().replace_residue(  6, new_ala2, False )
    
    pose.conformation().declare_chemical_bond(  6, "C", 7, "NE" )

    nres = pose.size()

    add_orn_cst( pose, nres, 1 )
    add_orn_cst( pose, 6, 7 )

    for ii in xrange(pose.size()):
        pose.set_phi(ii + 1, -150)
        pose.set_psi(ii + 1, 150)
        pose.set_omega(ii + 1, 180.0)

    # create move map for minimization
    mm = core.kinematics.MoveMap()
    mm.set_bb( True )
    mm.set_chi( True )
    mm.set_jump( 1, True )

    # create minimization mover
    minM = protocols.simple_moves.MinMover(mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, True)

    ####
    #Initial MC
    ####
    pert_sequence = protocols.moves.SequenceMover()
    pert_mc = protocols.moves.MonteCarlo(pose, score_fxn, 0.2)

    pert_pep_mm = core.kinematics.MoveMap()

    for i in xrange(pose.size()):
        pert_pep_mm.set_bb(i + 1)

    pert_pep_small = protocols.simple_moves.SmallMover(pert_pep_mm, 0.2, 1)
    pert_pep_small.angle_max('H', 2.0)
    pert_pep_small.angle_max('L', 2.0)
    pert_pep_small.angle_max('E', 2.0)
    pert_pep_shear = protocols.simple_moves.ShearMover(pert_pep_mm, 0.2, 1)
    pert_pep_shear.angle_max('H', 2.0)
    pert_pep_shear.angle_max('L', 2.0)
    pert_pep_shear.angle_max('E', 2.0)

    pert_pep_random = protocols.moves.RandomMover()
    pert_pep_random.add_mover(pert_pep_small, 1)
    pert_pep_random.add_mover(pert_pep_shear, 1)
    pert_pep_repeat = protocols.moves.RepeatMover(pert_pep_random, 100)

    pert_sequence.add_mover(pert_pep_repeat)
    pert_trial = protocols.moves.TrialMover(pert_sequence, pert_mc)

    # core.Size hbs_position = 1
    wt = 0.01
    for j in xrange(100):

        score_fxn.set_weight(core.scoring.atom_pair_constraint, wt)
        score_fxn.set_weight(core.scoring.angle_constraint, wt)
        score_fxn.set_weight(core.scoring.dihedral_constraint, wt)

        minM.apply(pose)

        pert_trial.apply(pose)

        #pert_mc.recover_low(pose)
        #pert_mc.show_state()
        pert_mc.boltzmann(pose)
        #TR<< "post mc.boltzmann" << std.endl
        #pert_mc.show_state()

        wt += 0.01
        print "After minimization, score is", (score_fxn)(pose) 

    pert_mc.recover_low(pose)
    minM.apply(pose)
    print "After minimization, score is", (score_fxn)(pose) 
    pose.dump_pdb( "out.pdb" )
