#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   lactamize.py
## @brief  Make any available lactam bridges in a PDB
## @author Andrew Watkins

def find_lactams(pose):
    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    restype_set = chm.residue_type_set('fa_standard').get_self_ptr()#.get_self_weak_ptr()

    resid_2_NC = {}
    for ii in xrange(pose.size()):
        if pose.residue_type(ii + 1).name3() == "LYS" or pose.residue_type(ii + 1).name3() == "ASP" or pose.residue_type(ii + 1).name3() == "GLU":
	    if pose.residue_type(ii + 1).has_property("BRANCH_POINT"): continue
            num_NC = len(resid_2_NC.keys())
	    resid_2_NC[ii + 1] = num_NC + 1
    
    if len(resid_2_NC.keys()) == 0: return
    print "Looking for lactam partners for", len(resid_2_NC.keys()), "residues."

    NC_2_resid = {}
    for ii in xrange(pose.size()):
        if ii + 1 in resid_2_NC:
            NC_2_resid[resid_2_NC[ii + 1]] = ii + 1


    # Create point graph of nbr_atoms of incomplete residues only.
    # Also, calculate maximum distance that a connected residue could be.
    pg = PointGraph()
    pg.set_num_vertices(len(resid_2_NC.keys()))
    maxrad = 0.0
    maxd = 1.33
    for ii in xrange(len(resid_2_NC.keys())):
        pg.get_vertex(ii + 1).data().xyz() = pose.residue(NC_2_resid[ii + 1]).atoms()[ii_res.nbr_atom()].xyz()
        if pose.residue(NC_2_resid[ii + 1]).nbr_radius() > maxrad: maxrad = pose.residue(NC_2_resid[ii + 1]).nbr_radius()

    # two Angstrom extra radius for finding bonds... very generous
    maxd += 2.0
    neighbor_cutoff = maxrad + maxd
    
    #find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff )

    # Iterate across neighbors of incomplete residues compare incomplete connection points against each other.
    for ii, ii_resid in NC_2_resid.iteritems():
        ii_res = pose.residue(ii_resid)
        if ii_res.has_property("BRANCH_POINT"): continue

        if ii_res.name().contains("N-conjugated"): continue
        if ii_res.name().contains("sidechain_carb"): continue

        jjatom = 0
        if ii_res.name3() == "LYS":
	    jjatom = ii_res.atom_index( "NZ" )
        elif ii_res.name3() == "ASP":
	    jjatom = ii_res.atom_index( "CG" )
        else:
	    jjatom = ii_res.atom_index( "CD" )

	best_match = 0.0
	best_match_resid = 0

	for ii_iter in pg.get_vertex(ii) in pg->get_vertex( ii ).upper_edge_list():

            neighb_id = ii_iter.upper_vertex()
	    # Disallow i, i+1!
	    neighb_resid = NC_2_resid[ neighb_id ]
            if neighb_resid == ii_resid + 1 or neighb_resid == ii_resid - 1: continue

	    neighb = pose.residue( neighb_resid )

            if neighb.has_property( "BRANCH_POINT" ): continue
            if neighb.name().contains("N-conjugated"): continue
            if neighb.name().contains("sidechain_carb"): continue

	    if ( ii_res.type().name3() == neighb.type().name3() ) continue
	    if ( ( ii_res.type().name3() == "GLU" && neighb.type().name3() == "ASP" ) ||
					( ii_res.type().name3() == "ASP" && neighb.type().name3() == "GLU" ) ) continue

			Size kkatom = 0
			if ( neighb.type().name3() == "LYS" ) {
				kkatom = neighb.type().atom_index( "NZ" )
			} else if ( neighb.type().name3() == "ASP" ) {
				kkatom = neighb.type().atom_index( "CG" )
			} else if ( neighb.type().name3() == "GLU" ) {
				kkatom = neighb.type().atom_index( "CD" )
			}

			// Calculate distances between expected atoms and actual atoms for these two connections.
			Distance kk_distance = ii_res.atom( jjatom ).xyz().distance( neighb.atom( kkatom ).xyz() )

			if ( best_match_resid == 0 || best_match > kk_distance ) {
				best_match = kk_distance
				best_match_resid = neighb_resid
			}
		}

		if best_match_resid == 0 ) {
			TR << "Failed to find a lactam partner for residue " << ii_resid << std::endl
                else:
			using namespace core
			using namespace id
			using namespace core::scoring
			using namespace core::scoring::constraints
			using namespace core::scoring::func

			TR << "Connecting residues: " << ii_resid << " ( " << ii_res.name()
			TR << " ) and " << best_match_resid << " ( " << pose.residue( best_match_resid ).name() << " )" << std::endl

			Size lys_id = 0
			Size glu_id = 0
			if ( ii_res.type().name3() == "LYS" ) {
				TR << "Replacing " << ii_resid << " with LYS:N-conjugated" << std::endl

				lys_id = ii_resid
				glu_id = best_match_resid

			} else {

				glu_id = ii_resid
				lys_id = best_match_resid

			}

			TR << "Replacing " << lys_id << " with LYS:N-conjugated" << std::endl

			replace_pose_residue_copying_existing_coordinates( pose, lys_id, residue_set_cap->name_map( "LYS:N-conjugated") )
			ResidueOP nr = ResidueFactory::create_residue( residue_set_cap->name_map("LYS:N-conjugated"), pose.residue( lys_id ), pose.conformation() )
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( lys_id ), *nr, pose.conformation() )
			nr->residue_connection_partner( pose.residue( lys_id ).n_polymeric_residue_connections() + 1, glu_id, pose.residue( glu_id ).n_polymeric_residue_connections() + 1 )
			pose.conformation().replace_residue( lys_id, *nr, false )

			std::string new_name = pose.residue(glu_id).name3() == "GLU"? "GLU:sidechain_carboxyl_conjugated" : "ASP:sidechain_carboxyl_conjugated"

			TR << "Replacing " << glu_id << " with "<<new_name << std::endl

			replace_pose_residue_copying_existing_coordinates( pose, glu_id, residue_set_cap->name_map( new_name ) )
			ResidueOP nr2 = ResidueFactory::create_residue( residue_set_cap->name_map(new_name), pose.residue( glu_id ), pose.conformation() )
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( glu_id ), *nr2, pose.conformation() )
			nr2->residue_connection_partner( pose.residue( glu_id ).n_polymeric_residue_connections() + 1, lys_id, pose.residue( lys_id ).n_polymeric_residue_connections() + 1 )
			pose.conformation().replace_residue( glu_id, *nr2, false )



			// add constraints
			HarmonicFuncOP harm( new core::scoring::func::HarmonicFunc( 1.33, 0.05 ) )
			CircularHarmonicFuncOP ang180( new core::scoring::func::CircularHarmonicFunc( 3.14159, 0.05 ) )
			CircularHarmonicFuncOP ang120( new core::scoring::func::CircularHarmonicFunc( 3.14159*2.0/3.0, 0.05 ) )
			AtomID const lysHZ( pose.residue( lys_id ).atom_index( "1HZ" ), best_match_resid )
			AtomID const lysNZ( pose.residue( lys_id ).atom_index( "NZ" ), best_match_resid )
			AtomID const lysCD( pose.residue( lys_id ).atom_index( "CD" ), best_match_resid )
			if ( pose.residue( glu_id ).name3() == "GLU" ) {
				AtomID const CD( pose.residue( glu_id ).atom_index( "CD"  ), ii_resid )
				AtomID const OE( pose.residue( glu_id ).atom_index( "OE1"  ), ii_resid )
				AtomID const CG( pose.residue( glu_id ).atom_index( "CG"  ), ii_resid )
				pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( lysHZ, lysNZ, CD, OE, ang180 ) ) )
				pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( lysNZ, CD, OE, CG, ang180 ) ) )
				pose.add_constraint( AngleConstraintOP( new AngleConstraint( lysNZ, CD, OE, ang120 ) ) )
				pose.add_constraint( AngleConstraintOP( new AngleConstraint( lysHZ, lysNZ, CD, ang120 ) ) )

				pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint( lysNZ, CD, harm ) ) )
				pose.conformation().declare_chemical_bond(  glu_id, "CD", lys_id, "NZ" )

			} else {
				AtomID const CG( pose.residue( glu_id ).atom_index( "CG"  ), ii_resid )
				AtomID const OD( pose.residue( glu_id ).atom_index( "OD1"  ), ii_resid )
				AtomID const CB( pose.residue( glu_id ).atom_index( "CB"  ), ii_resid )
				pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( lysHZ, lysNZ, CG, OD, ang180 ) ) )
				pose.add_constraint( DihedralConstraintOP( new DihedralConstraint( lysNZ, CG, OD, CB, ang180 ) ) )
				pose.add_constraint( AngleConstraintOP( new AngleConstraint( lysNZ, CG, OD, ang120 ) ) )
				pose.add_constraint( AngleConstraintOP( new AngleConstraint( lysHZ, lysNZ, CG, ang120 ) ) )

				pose.add_constraint( AtomPairConstraintOP( new AtomPairConstraint( lysNZ, CG, harm ) ) )
				pose.conformation().declare_chemical_bond(  glu_id, "CG", lys_id, "NZ" )

			}

		}

	}
}

// find all ostensible lactams with bad bond distances and revert those residues...
void revert_lactams( Pose & pose ) {

	core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )

	for ( Size ii = 1 ii <= pose.size() ++ii ) {
		if ( pose.residue(ii).type().name().find("N-conjugated")!= std::string::npos || pose.residue(ii).type().name().find("sidechain_carb")!= std::string::npos ) {

			Size other_res = pose.residue( ii ).residue_connection_partner( pose.residue( ii ).n_polymeric_residue_connections() + 1 )
			TR << "Testing reversion of " << ii << " and  " << other_res << std::endl
			TR << "Measuring distance from " << pose.residue(ii).atom_name( pose.residue(ii).residue_connect_atom_index( pose.residue( ii ).n_polymeric_residue_connections() + 1 ) ) << " to " << pose.residue(other_res).atom_name( pose.residue(other_res).residue_connect_atom_index( pose.residue( other_res).n_polymeric_residue_connections() + 1 ) ) << std::endl
			Distance dist = pose.residue(ii).atom( pose.residue(ii).residue_connect_atom_index( pose.residue( ii ).n_polymeric_residue_connections() + 1 ) ).xyz().distance(
				pose.residue(other_res).atom( pose.residue(other_res).residue_connect_atom_index( pose.residue( other_res).n_polymeric_residue_connections() + 1 ) ).xyz() )
			TR << "Distance is " << dist << std::endl
			if ( dist > 1.6 ) {
				// revert both!

				Size lys_id = 0
				Size glu_id = 0
				if ( pose.residue(ii).type().name3() == "LYS" ) {
					lys_id = ii
					glu_id = other_res
				} else {
					glu_id = ii
					lys_id = other_res
				}


				core::scoring::constraints::ConstraintCOPs cs = pose.constraint_set()->get_all_constraints()
				for ( Size i = 1 i <= cs.size() i++ ) {
					Constraint const & other_cst = *cs[i]
					if ( !dynamic_cast< AtomPairConstraint const * > ( &other_cst ) ) {
						continue
					}

					auto const & constraint_i( static_cast< AtomPairConstraint const & > (other_cst) )

					if ( constraint_i.atom(1) == AtomID(pose.residue(glu_id).atom_index("NZ"), glu_id) && constraint_i.atom(2) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CD" : "CG" ), glu_id) ) {
						pose.remove_constraint( cs[i], true )
					}
				}
				for ( Size i = 1 i <= cs.size() i++ ) {
					Constraint const & other_cst = *cs[i]
					if ( !dynamic_cast< AngleConstraint const * > ( &other_cst ) ) {
						continue
					}

					auto const & constraint_i( static_cast< AngleConstraint const & > (other_cst) )

					if ( constraint_i.atom(1) == AtomID(pose.residue(lys_id).atom_index( "NZ" ), lys_id) && constraint_i.atom(2) == AtomID( pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CD" : "CG" ), glu_id) && constraint_i.atom(3) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "OE1" : "OD1" ), glu_id) ) {
						pose.remove_constraint( cs[i], true )
					}
					if ( constraint_i.atom(1) == AtomID(pose.residue(lys_id).atom_index( "1HZ" ), lys_id) && constraint_i.atom(2) == AtomID(pose.residue(lys_id).atom_index( "NZ" ),lys_id) && constraint_i.atom(3) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CD" : "CG" ),glu_id) ) {
						pose.remove_constraint( cs[i], true )
					}
				}
				for ( Size i = 1 i <= cs.size() i++ ) {
					Constraint const & other_cst = *cs[i]
					if ( !dynamic_cast< DihedralConstraint const * > ( &other_cst ) ) {
						continue
					}

					auto const & constraint_i( static_cast< DihedralConstraint const & > (other_cst) )

					if ( constraint_i.atom(1) == AtomID(pose.residue(lys_id).atom_index( "1HZ" ),lys_id) && constraint_i.atom(2) == AtomID(pose.residue(lys_id).atom_index( "NZ" ),lys_id) && constraint_i.atom(3) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CD" : "CG" ),glu_id)  &&  constraint_i.atom(4) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "OE1" : "OD1" ),glu_id) ) {
						pose.remove_constraint( cs[i], true )
					}
					if ( constraint_i.atom(1) == AtomID(pose.residue(lys_id).atom_index( "NZ" ),lys_id) &&  constraint_i.atom(2) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CD" : "CG" ),glu_id) &&  constraint_i.atom(3) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "OE1" : "OD1" ),glu_id) &&  constraint_i.atom(4) == AtomID(pose.residue(glu_id).atom_index( pose.residue(glu_id).name3() == "GLU"? "CG" : "CB" ),glu_id)  ) {
						pose.remove_constraint( cs[i], true )
					}
				}


				TR << "Reverting " << lys_id << " with LYS" << std::endl

				replace_pose_residue_copying_existing_coordinates( pose, lys_id, residue_set_cap->name_map( "LYS") )
				ResidueOP nr = ResidueFactory::create_residue( residue_set_cap->name_map("LYS"), pose.residue( lys_id ), pose.conformation() )
				core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( lys_id ), *nr, pose.conformation() )
				nr->residue_connection_partner( pose.residue( lys_id ).n_polymeric_residue_connections() + 1, glu_id, pose.residue( glu_id ).n_polymeric_residue_connections() + 1 )
				pose.conformation().replace_residue( lys_id, *nr, false )

				std::string new_name = pose.residue(glu_id).name3() == "GLU"? "GLU" : "ASP"

				TR << "Reverting " << glu_id << " with "<<new_name << std::endl

				replace_pose_residue_copying_existing_coordinates( pose, glu_id, residue_set_cap->name_map( new_name ) )
				ResidueOP nr2 = ResidueFactory::create_residue( residue_set_cap->name_map(new_name), pose.residue( glu_id ), pose.conformation() )
				core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( glu_id ), *nr2, pose.conformation() )
				nr2->residue_connection_partner( pose.residue( glu_id ).n_polymeric_residue_connections() + 1, lys_id, pose.residue( lys_id ).n_polymeric_residue_connections() + 1 )
				pose.conformation().replace_residue( glu_id, *nr2, false )

			}
		}
	}
}

int
main( int argc, char* argv[] )
{
	try {
		devel::init(argc, argv)

		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function()

		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1 )
		score_fxn->set_weight( core::scoring::angle_constraint, 1.0 )
		score_fxn->set_weight( core::scoring::dihedral_constraint, 1 )//10.0 )

		Pose pose
		std::string filename = option[in::file::s].value()[1]

		import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file)
		TR << "Importing pose from " << filename << std::endl
		filename.erase( filename.size()-4 )

		// pack rotamers loop
		for ( Size ii = 1 ii <= 100 ++ii ) {
			std::string outfile = filename
			outfile += "_" + utility::to_string( ii ) + ".pdb"


			pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( pose ))
			base_packer_task->set_bump_check( false )
			base_packer_task->initialize_from_command_line()
			base_packer_task->or_include_current( true )
			for ( Size ii = 1 ii <= pose.size() ++ii ) {
				if ( pose.residue(ii).type().name().find("N-conjugated")!= std::string::npos || pose.residue(ii).type().name().find("sidechain_carb")!= std::string::npos ) {
					base_packer_task->nonconst_residue_task(ii).prevent_repacking()
				} else {
					base_packer_task->nonconst_residue_task(ii).restrict_to_repacking()
				}
			}
			pack::pack_rotamers( pose, *score_fxn, base_packer_task )


			kinematics::MoveMapOP pert_mm( new kinematics::MoveMap() )
			for ( Size ii = 1 ii <= pose.size() ++ii ) {
				if ( pose.residue(ii).type().name().find("N-conjugated")!= std::string::npos || pose.residue(ii).type().name().find("sidechain_carb")!= std::string::npos ) {
					pert_mm->set_chi( ii, true )
				} else if ( ii > 1 && (pose.residue(ii-1).type().name().find("N-conjugated")!= std::string::npos || pose.residue(ii-1).type().name().find("sidechain_carb")!= std::string::npos ) ) {
					pert_mm->set_chi( ii, true )
				} else if ( ii < pose.size() && ( pose.residue(ii+1).type().name().find("N-conjugated")!= std::string::npos || pose.residue(ii+1).type().name().find("sidechain_carb")!= std::string::npos ) ) {
					pert_mm->set_chi( ii, true )
				}
			}
			protocols::simple_moves::MinMoverOP min_mover( new simple_moves::MinMover( pert_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) )
			min_mover->cartesian( score_fxn->get_weight( cart_bonded )  > 0 )

			// neighbor graph
			find_lactams( pose )
			min_mover->apply( pose )
			revert_lactams( pose )

			pose.dump_pdb( outfile )
		}
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl
		return -1
	}
	return 0

}//main

