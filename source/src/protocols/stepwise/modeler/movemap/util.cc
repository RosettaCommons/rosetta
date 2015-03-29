// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/movemap/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
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

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.movemap.util" );

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
namespace modeler {
namespace movemap {

/////////////////////////////////////////////////////////////////////////////////////////////
// used by legacy StepWiseProteinMinimizer.cc
void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 core::pose::Pose const & pose,
														 utility::vector1< Size > const & working_minimize_res,
														 bool const move_takeoff_torsions /* = true */ ){
	toolbox::AllowInsertOP allow_insert( new toolbox::AllowInsert( pose ) );
	figure_out_stepwise_movemap( mm, allow_insert, pose, working_minimize_res, move_takeoff_torsions );
}

/////////////////////////////////////////////////////////////////////////////////////////////
// used by legacy StepWiseRNA_Minimizer.cc
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
	allow_insert = toolbox::AllowInsertOP( new toolbox::AllowInsert( pose ) );
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

			if ( check_sample_sugar_in_full_model_info( pose, n ) ) allow_insert->set_sugar( n, pose, true );
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
	//	output_movemap( mm, pose, TR );
}

/////////////////////////////////////////////////////////
// Useful checks...
// no longer appears in use -- deprecate in 2015
// if still not in use then.
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


} //movemap
} //modeler
} //stepwise
} //protocols
