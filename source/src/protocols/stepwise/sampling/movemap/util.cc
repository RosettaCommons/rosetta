// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/movemap/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/movemap/util.hh>
#include <protocols/stepwise/sampling/output_util.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.movemap/util" );

/////////////////////////////////////////////////////////////////////////////////////////////
//
// All movemap setting for stepwise assembly/monte carlo is now tucked into here.
//  Might consider making this a class.
//  Also might move "core" figure_out_stepwise_movemap() routine that goes from
//    AllowInsert to MoveMap into AllowInsert. Its pretty general and good.
//
//   -- rhiju, 2014
//
/////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {

/////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 core::pose::Pose const & pose,
														 utility::vector1< Size > const & working_minimize_res,
														 bool const move_takeoff_torsions /* = true */ ){
	toolbox::AllowInsertOP allow_insert = new toolbox::AllowInsert( pose );
	figure_out_stepwise_movemap( mm, allow_insert, pose, working_minimize_res, move_takeoff_torsions );
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 toolbox::AllowInsertOP & allow_insert,
														 core::pose::Pose const & pose,
														 utility::vector1< Size > const & working_fixed_res,
														 utility::vector1< Size > const & working_extra_minimize_res,
														 bool const move_takeoff_torsions /* = true */ ){

	utility::vector1< core::Size > working_minimize_res;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( !working_fixed_res.has_value( n ) || working_extra_minimize_res.has_value( n ) ) working_minimize_res.push_back( n );
	}
	allow_insert = new toolbox::AllowInsert( pose );
	figure_out_stepwise_movemap( mm, allow_insert, pose, working_minimize_res, move_takeoff_torsions );
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 toolbox::AllowInsertOP allow_insert,
														 core::pose::Pose const & pose,
														 utility::vector1< Size > const & working_minimize_res,
														 bool const move_takeoff_torsions /* = true */ ) {
	using namespace core::id;
	Size const nres( pose.total_residue() );
	//	TR << "WORKING MINIMIZE RES " << working_minimize_res << std::endl;
	for ( Size n = 1; n <= nres; n++ ) {
		if ( working_minimize_res.has_value( n ) ) allow_insert->set( n, true );
		if ( !working_minimize_res.has_value( n ) ) {
			if ( allow_insert->get( n ) ) allow_insert->set( n, false );
		}
		if ( pose.residue_type( n ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
			allow_insert->set( NamedAtomID( "OVL1",  n), pose, true );
			allow_insert->set( NamedAtomID( "OVL2",  n), pose, true );
		}
		if ( pose.residue_type( n ).has_variant_type( chemical::CUTPOINT_UPPER ) ){
			allow_insert->set( NamedAtomID( "OVU1", n ), pose, true );
		}
		if ( pose.residue_type( n ).is_RNA() ){
			// there's a problem with looking at VIRTUAL PHOSPHATE -- it can change if 5' packing phosphate is 'packed' in.
			// instead let it move -- but rely on the fact that there are no score terms computed (so derivative will be zero)
			//			if ( pose.residue_type(n).has_variant_type( chemical::VIRTUAL_PHOSPHATE ) ) allow_insert->set_phosphate( n, pose, false );

			// following are 'conventions' that need to be hard-coded -- for RNA want
			// entire suites (not just connection torsions ) to move between residues in different domains
			if ( pose.residue_type(n).has_variant_type( chemical::CUTPOINT_UPPER    ) ) allow_insert->set_phosphate( n, pose, true );

			if ( n > 1 ) {
				if ( move_takeoff_torsions && pose.residue_type( n-1 ).is_RNA() ){
					if ( working_minimize_res.has_value( n-1 ) )	allow_insert->set_phosphate( n, pose, true );
					if ( allow_insert->get_domain( NamedAtomID( " O3'", n-1 ), pose ) !=
							 allow_insert->get_domain( NamedAtomID( " P  ", n )  , pose ) ) allow_insert->set_phosphate( n, pose, true );
				}
			}
		}
		if ( pose.residue_type(n).is_protein() ){
			if ( move_takeoff_torsions ){
				if ( n > 1 && pose.residue_type( n-1 ).is_protein() && working_minimize_res.has_value( n-1 ) &&
						 pose.residue_type( n ).has( " H  " ) /* not there for Pro*/ ) {
					allow_insert->set( NamedAtomID(" H  ", n ), pose, true );
				}
				if ( n < nres && pose.residue_type( n+1 ).is_protein() && working_minimize_res.has_value( n+1 ) ) {
					allow_insert->set( NamedAtomID(" O  ", n ), pose, true );
				}
			} else { // argh, proteins. This is to match old KIC runs. Perhaps should just deprecate.
				if ( n > 1 && pose.residue_type( n-1 ).is_protein() &&
						 !working_minimize_res.has_value( n-1 ) ) {
					allow_insert->set_domain( NamedAtomID(" N  ", n ), pose,
																		allow_insert->get_domain( NamedAtomID( " C  ", n-1 ), pose ) );
					allow_insert->set_domain( NamedAtomID(" CA ", n ), pose,
																		allow_insert->get_domain( NamedAtomID( " C  ", n-1 ), pose ) );
				}
				if ( n < nres && pose.residue_type( n+1 ).is_protein() &&
						 !working_minimize_res.has_value( n+1 ) ) {
					allow_insert->set_domain( NamedAtomID(" C  ", n ), pose,
																		allow_insert->get_domain( NamedAtomID( " N  ", n+1 ), pose ) );
				}
			}
		}
	}

	allow_insert->setup_movemap( mm, pose );
	//output_movemap( mm, pose, TR );
}

////////////////////////////////////////////////////////////////////////////////////
// used in protein SWA -- should deprecate soon in favor of general movemap setup.
////////////////////////////////////////////////////////////////////////////////////
void
figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose,
														utility::vector1< core::Size > const & fixed_res,
														bool const move_takeoff_torsions,
														bool const move_jumps_between_chains )
{

	using namespace core::id;

	Size const nres( pose.total_residue() );

	ObjexxFCL::FArray1D< bool > allow_insert( nres, true );
	for (Size i = 1; i <= fixed_res.size(); i++ ) allow_insert( fixed_res[ i ] ) = false;
	for (Size n = 1; n <= pose.total_residue(); n++ ){
		if ( pose.residue( n ).has_variant_type( "VIRTUAL_RESIDUE" ) ) allow_insert( n ) = false;
	}

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	for  (Size i = 1; i <= nres; i++ )  {

		if ( !move_takeoff_torsions && !allow_insert(i) ) continue; // don't allow, e.g., psi/omega of residue before loop to move.

		utility::vector1< TorsionID > torsion_ids;

		for ( Size torsion_number = 1; torsion_number <= pose.residue( i ).mainchain_torsions().size(); torsion_number++ ) {
			torsion_ids.push_back( TorsionID( i, id::BB, torsion_number ) );
		}
		for ( Size torsion_number = 1; torsion_number <= pose.residue( i ).nchi(); torsion_number++ ) {
			torsion_ids.push_back( TorsionID( i, id::CHI, torsion_number ) );
		}

		for ( Size n = 1; n <= torsion_ids.size(); n++ ) {

			TorsionID const & torsion_id  = torsion_ids[ n ];

			id::AtomID id1,id2,id3,id4;
			bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
			if (fail) continue;

			// If there's any atom that is in a moving residue by this torsion, let the torsion move.
			//  should we handle a special case for cutpoint atoms? I kind of want those all to move.
			if ( !allow_insert( id1.rsd() ) && !allow_insert( id2.rsd() ) && !allow_insert( id3.rsd() )  && !allow_insert( id4.rsd() ) ) continue;
			mm.set(  torsion_id, true );

		}
	}

	utility::vector1< Size > chain_index;
	Size chain_number( 0 );
	for (Size n = 1; n <= pose.total_residue(); n++ ){
		if ( pose.residue_type( n ).is_lower_terminus()  && !pose.residue_type( n ).has_variant_type( chemical::N_ACETYLATION ) ) chain_number++;
		chain_index.push_back( chain_number );
	}

	for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

		if ( allow_insert( jump_pos1 ) || allow_insert( jump_pos2 ) ) 	 {
			mm.set_jump( n, true );
			//				TR  << "allow_insert ALLOWING JUMP " << n << " to move. It connects " << jump_pos1 << " and " << jump_pos2 << "." << std::endl;
		}

		if ( move_jumps_between_chains ){
			if ( chain_index[ jump_pos1 ] != chain_index[ jump_pos2 ] ){
				utility_exit_with_message( "is this move_jumps_between_chains_ option still in use? The mm.set_jump() command is gone." );
				TR << "move_jumps_between_chains ALLOWING JUMP " << n << " to move. It connects " << jump_pos1 << " and " << jump_pos2 << "." << std::endl;
			}
		}

	}

}


/////////////////////////////////////////////////////////
// Useful checks...
/////////////////////////////////////////////////////////
void
check_move_map_against_working_parameters( core::pose::Pose const & pose ,
																					 core::kinematics::MoveMapCOP minimize_move_map,
																					 working_parameters::StepWiseWorkingParametersCOP working_parameters ){
	using namespace core::id;

	utility::vector1< Size > suites_that_must_be_minimized;
	////////////////////////////////////////////////////////////////////////////////////
	// make sure moving dof is actually being minimized.
	Size const rebuild_suite = working_parameters->working_moving_suite(); // is this set up correctly?
	if ( rebuild_suite > 0 && !working_parameters->floating_base() ) {
		suites_that_must_be_minimized.push_back( rebuild_suite );
	}
	utility::vector1< Size > cutpoint_closed = figure_out_moving_cutpoints_closed( pose, working_parameters->working_moving_partition_res() );
	for ( Size i = 1; i <= cutpoint_closed.size(); i++ ) suites_that_must_be_minimized.push_back( cutpoint_closed[i] );

	// last, but not least, there might be some information in the domain map. Note
	// that generally we could instead replace working_fixed_res with an inputted domain map.
	// that is, get rid of working_fixed_res & minimize_res and instead have a local fixed_domain_map,
	// which can instead be updated by set_working_fixed_res.
	// how to tell Modeler to *not* minimize additional suites?
	suites_that_must_be_minimized = merge_vectors( suites_that_must_be_minimized, get_domain_boundary_suites( pose ) );

	TR.Debug << "SUITES_THAT_MUST_BE_MINIMIZED " << suites_that_must_be_minimized << std::endl;
	for ( Size n = 1; n <= suites_that_must_be_minimized.size(); n++ ){
		Size const suite_num = suites_that_must_be_minimized[ n ];
		if ( pose.residue_type( n ).is_RNA() ) {
			runtime_assert( minimize_move_map->get( TorsionID( suite_num,   id::BB, 5 ) ) ); // epsilon
			runtime_assert( minimize_move_map->get( TorsionID( suite_num,   id::BB, 6 ) ) ); // zeta
			runtime_assert( minimize_move_map->get( TorsionID( suite_num+1, id::BB, 1 ) ) ); // alpha
			runtime_assert( minimize_move_map->get( TorsionID( suite_num+1, id::BB, 2 ) ) ); // beta
			runtime_assert( minimize_move_map->get( TorsionID( suite_num+1, id::BB, 3 ) ) ); // gamma
		} else if ( pose.residue_type( n ).is_protein() ) {
			runtime_assert( minimize_move_map->get( TorsionID( suite_num,   id::BB, 2 ) ) ); // epsilon
			runtime_assert( minimize_move_map->get( TorsionID( suite_num,   id::BB, 3 ) ) ); // zeta
			runtime_assert( minimize_move_map->get( TorsionID( suite_num+1, id::BB, 1 ) ) ); // alpha
		}
	}

	if ( working_parameters->floating_base() ) {
		Size const rebuild_res = working_parameters->working_moving_res();
		runtime_assert( minimize_move_map->get_jump( pose.fold_tree().get_jump_that_builds_residue( rebuild_res ) ) );
	}
}


} //sampling
} //stepwise
} //protocols
