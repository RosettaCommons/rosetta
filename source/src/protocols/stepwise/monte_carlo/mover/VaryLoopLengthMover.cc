// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.cc
/// @brief In stepwise design, vary desired loop lengths by updating FullModelParameters
/// @author Rhiju Das (rhiju@stanford.edu)

#include <protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/util.hh>

#include <protocols/stepwise/modeler/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.VaryLoopLengthMover" );


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

VaryLoopLengthMover::VaryLoopLengthMover():
	protocols::moves::Mover( "VaryLoopLengthMover" )
{
}

VaryLoopLengthMover::~VaryLoopLengthMover(){}

VaryLoopLengthMover::VaryLoopLengthMover( VaryLoopLengthMover const & src ):
	protocols::moves::Mover( src )
{

}

///////////////////////////////////////////////////////////////////////////
std::string
VaryLoopLengthMover::get_name() const {
	return "VaryLoopLengthMover";
}

///////////////////////////////////////////////////////////////////////////
void
VaryLoopLengthMover::apply( core::pose::Pose &  ){
	// no op -- must supply swa_move.
}

///////////////////////////////////////////////////////////////////////////
//
// Example -- suppose we have a loop that can be up to 5 residues long, but for
//  now we're only assuming it has 3 residues.
//
// The parent_full_model_parameters will know about the 5-residue loop:
//
//  1  2  3  4  5  6  7  8  9  10
//        L  L  L  L  L
//
// The full_model_parameters will simulate having the 3-residue loop version:
//
//  1  2  3  4  5  8  9  10   <-- slice_res_list, i.e. what each full_model residue is called in parent_full_model_numbering
//        L  L  L
//
// This class will update slice_res_list & full_model_parameters if we decrease or
//   increase desired loop length. E.g. upon  ( 3 DELETE_LOOP_RES BOND_TO_PREVIOUS 2 ), should get:
//
//  1  2  3  4  8  9  10   <-- slice_res_list
//        L  L
//
// Need to also propagate change into res_list held by each pose.
//
void
VaryLoopLengthMover::apply( core::pose::Pose & pose, StepWiseMove const & swa_move ) {
	using namespace core::pose::full_model_info;
	FullModelParametersCOP full_model_parameters( const_full_model_info( pose ).full_model_parameters() );
	FullModelParametersCOP parent_full_model_parameters( full_model_parameters->parent_full_model_parameters() );
	runtime_assert( parent_full_model_parameters != 0 );

	// which parent_full_model residues are being modeled in current full_model_parameters:
	utility::vector1< Size > const & slice_res_list( full_model_parameters->slice_res_list() );

	if ( swa_move.move_type() == DELETE_LOOP_RES ) {
		Size const & delete_res = swa_move.moving_res();
		//  runtime_assert( !all_res_list.has_value( delete_res ) );
		Size const delete_res_in_parent_full_model_numbering = slice_res_list[ delete_res ];

		utility::vector1< Size > slice_res_list_new;
		for ( Size n = 1; n <= slice_res_list.size(); n++ ) {
			if ( slice_res_list[ n ] != delete_res_in_parent_full_model_numbering ) {
				slice_res_list_new.push_back( slice_res_list[ n ] ) ;
			}
		}
		runtime_assert( slice_res_list_new.size() == slice_res_list.size() - 1 );

		FullModelParametersCOP full_model_parameters_new = parent_full_model_parameters->slice( slice_res_list_new );
		update_full_model_parameters( pose, full_model_parameters_new );
	} else {
		runtime_assert( swa_move.move_type() == ADD_LOOP_RES );
		if ( swa_move.attachment_type() == BOND_TO_PREVIOUS ) {
			runtime_assert( swa_move.moving_res() == swa_move.attached_res() + 1 );
		} else {
			runtime_assert( swa_move.attachment_type() == BOND_TO_NEXT );
			runtime_assert( swa_move.moving_res() == swa_move.attached_res() - 1 );
		}
		Size const attached_res_in_parent_full_model_numbering = slice_res_list[ swa_move.attached_res() ];
		Size const add_res_in_parent_full_model_numbering = swa_move.attachment_type() == BOND_TO_PREVIOUS ?
			attached_res_in_parent_full_model_numbering + 1 :
			attached_res_in_parent_full_model_numbering - 1;
		runtime_assert( !slice_res_list.has_value( add_res_in_parent_full_model_numbering ) );

		utility::vector1< Size > slice_res_list_new;
		for ( Size n = 1; n <= slice_res_list.size(); n++ ) {
			slice_res_list_new.push_back( slice_res_list[ n ] ) ;
			if ( slice_res_list[ n ] == add_res_in_parent_full_model_numbering - 1 ) {
				slice_res_list_new.push_back( add_res_in_parent_full_model_numbering );
			}
		}
		runtime_assert( slice_res_list_new.size() == slice_res_list.size() + 1 );

		FullModelParametersCOP full_model_parameters_new = parent_full_model_parameters->slice( slice_res_list_new );
		update_full_model_parameters( pose, full_model_parameters_new );
	}

	modeler::fix_up_jump_atoms_and_residue_type_variants( pose );
}

// this should recurse down through other_poses, updating their full_model_parameters & slice_res_list.
// ... and it needs to update res_list in each pose.
void
VaryLoopLengthMover::update_full_model_parameters( pose::Pose & pose, pose::full_model_info::FullModelParametersCOP full_model_parameters_new ) const
{
	using namespace core::pose::full_model_info;
	using namespace utility;
	FullModelParametersCOP full_model_parameters( const_full_model_info( pose ).full_model_parameters() );

	// better be derived from same parent full_model with maximal loop lengths.
	runtime_assert( full_model_parameters->parent_full_model_parameters() != 0 );
	runtime_assert( full_model_parameters->parent_full_model_parameters() ==
		full_model_parameters_new->parent_full_model_parameters() );

	vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
	vector1< Size > const & slice_res_list = full_model_parameters->slice_res_list();
	vector1< Size > const & slice_res_list_new = full_model_parameters_new->slice_res_list();
	vector1< Size > res_list_new;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		// note that all instantiated residues in the pose  need to be carried over --
		// that is, DELETE_LOOP_RES and ADD_LOOP_RES moves cannot delete or add instantiated loop residues,
		// just potential loop residues....
		Size const res_in_parent_full_model_numbering( slice_res_list[ res_list[ n ] ] );
		runtime_assert( slice_res_list_new.has_value( res_in_parent_full_model_numbering  ) );
		res_list_new.push_back( slice_res_list_new.index( res_in_parent_full_model_numbering ) );
	}

	FullModelInfoOP full_model_info_new = const_full_model_info( pose ).clone_info();
	full_model_info_new->set_res_list( res_list_new );
	full_model_info_new->set_full_model_parameters( full_model_parameters_new );
	set_full_model_info( pose, full_model_info_new );
	runtime_assert( check_full_model_info_OK( pose ) );
	update_pose_objects_from_full_model_info( pose );

	// recurse through other poses.
	vector1< pose::PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ) {
		update_full_model_parameters( *other_pose_list[ n ], full_model_parameters_new );
	}
}


} //protocols
} //stepwise
} //monte_carlo
} //mover


