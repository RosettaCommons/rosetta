#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   create_a3b_hbs.py
## @brief  Scan the sidechain angles of a beta-AA
## @author Andrew Watkins

def set_up_movemap(pose):
    """
    Set up a movemap that allows every DOF to minimize except for 
    certain chis.
    """
    mm = core.kinematics.MoveMap()
    for ii in xrange(pose.residue_type(1).natoms()):
        jj = ii + 1
        mm.set(core.id.DOF_ID(core.id.AtomID(jj, 1), core.id.D), True)
        mm.set(core.id.DOF_ID(core.id.AtomID(jj, 1), core.id.THETA), True)

    for ii in xrange(pose.residue_type(1).nchi()):
        jj = ii + 1
        if pose.residue_type(1).is_proton_chi(jj):
            mm.set(core.id.TorsionID( 1, core.id.CHI, 1), True)


    mm.show()
    return mm
    #mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "CO" ), 1 ), PHI ), false );
    #mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "C" ), 1 ), PHI ), false );
    #mm->set( DOF_ID( AtomID( pose.residue_type( 1 ).atom_index( "NM" ), 1 ), PHI ), false );

def scan_chi(pose, score_fxn, mm, extension):
    minmover = protocols.simple_moves.MinMover(mm, score_fxn, "linmin_iterated", 0.001, True)
    
    #utility::vector1< utility::vector1< Real > > results;
    if pose.residue_type( 1 ).nchi() == 1:
        for chi1 in xrange(0, 360, 10):
	    #utility::vector1< Real > row;
            pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.CHI, 1 ), chi1 )
            # score the pose
            #Real orig_ener( ( *score_fxn )( pose ) );
	    minmover.apply(pose)
	    pose.dump_pdb("pose_{chi1}_{ext}.pdb".format(chi1=chi1, ext=extension))
	    min_ener = pose.energies().total_energy()
	    #row.push_back(min_ener);
	    #results.push_back( row );
	#print results
    elif pose.residue_type( 1 ).nchi() == 2:
        for chi1 in xrange(0, 360, 10):
            for chi2 in xrange(0, 360, 10):
	        #utility::vector1< Real > row;
                pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.CHI, 1 ), chi1 )
                pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.CHI, 2 ), chi2 )
                # score the pose
                #Real orig_ener( ( *score_fxn )( pose ) );
	        minmover.apply(pose)
	        pose.dump_pdb("pose_{chi1}_{chi2}_{ext}.pdb".format(chi1=chi1, chi2=chi2, ext=extension))
	        min_ener = pose.energies().total_energy()
	        #row.push_back(min_ener);
	        #results.push_back( row );
	    #print results

if __name__ == '__main__':
    from pyrosetta import *
    from pyrosetta.rosetta import *
    import sys

    init()

    score_fxn = core.scoring.get_score_function()
    name3 = sys.argv[1]
    #assert( name3 in ['B3A', 'B3C', 'B3D', 'B3E', 'B3F', 'B3G', 'B3H', 'B3I', 'B3K', 'B3L', 'B3M', 'B3N', 'B3O', 'B3P', 'B3Q', 'B3R', 'B3S', 'B3T', 'B3V', 'B3W', 'B3Y'] )

    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    restype_set = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()
    seq = 'X[{name3}:AcetylatedNtermProteinFull:MethylatedCtermProteinFull]'.format(name3=name3) 
    pose = pyrosetta.pose_from_sequence(seq)
    mm = set_up_movemap(pose)
    
    # Extended
    pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 1 ), -140 );
    pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 2 ),   75 );
    #pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 3 ),   65 );

    scan_chi( pose, score_fxn, mm, "ext" );
    
    # Helical
    pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 1 ), -140 );
    pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 2 ),   60 );
    #pose.conformation().set_torsion( core.id.TorsionID( 1, core.id.BB, 3 ), -120 );

    scan_chi( pose, score_fxn, mm, "hel" );
