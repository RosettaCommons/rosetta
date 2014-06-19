// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>

#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <protocols/stepwise/sampling/util.hh>

#include <basic/Tracer.hh>

using namespace core::chemical;
using namespace core;

static basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.TransientCutpointHandler" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

	//Constructor
	TransientCutpointHandler::TransientCutpointHandler( Size const sample_res ):
		sample_suite_  ( sample_res - 1 ), //sampled nucleoside
		cutpoint_suite_( sample_res ), // suite at which to make cutpoint
		move_jump_points_away_( false )
	{}

	//Constructor
	TransientCutpointHandler::TransientCutpointHandler( Size const sample_suite, Size const cutpoint_suite ):
		sample_suite_( sample_suite ), //sampled nucleoside
		cutpoint_suite_( cutpoint_suite ), // suite at which to make cutpoint
		move_jump_points_away_( false )
	{}

	//Destructor
	TransientCutpointHandler::~TransientCutpointHandler()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	void
	TransientCutpointHandler::put_in_cutpoints( core::pose::Pose & viewer_pose ){

		Pose pose = viewer_pose; // prevent some conflicts with graphics. Note potential slowdown.

		fold_tree_save_ = pose.fold_tree();

		// create reasonable fold tree
		prepare_fold_tree_for_erraser( pose );

		pose::correctly_add_cutpoint_variants( pose, cutpoint_suite_ );

		viewer_pose = pose;
	}

	/////////////////////////////////////////////////////////////////////////
	void
	TransientCutpointHandler::prepare_fold_tree_for_erraser( core::pose::Pose & pose ){

		using namespace core::kinematics;
		using namespace core::chemical::rna;

		FoldTree f = pose.fold_tree();

		//update_fixed_res_and_minimize_res( pose );
		utility::vector1< Size > sample_res_list = minimize_res_;

		// figure out jump points that bracket the cut.
		Size jump_start( sample_suite_    );
		Size jump_end( cutpoint_suite_ + 1);

		if ( move_jump_points_away_ ){
			while ( jump_start > 1 && minimize_res_.has_value( jump_start ) && !f.is_cutpoint( jump_start - 1 ) )          jump_start--;
			while ( jump_end < pose.total_residue() && minimize_res_.has_value( jump_end ) && !f.is_cutpoint( jump_end ) ) jump_end++;
		}

		Size const cutpoint = cutpoint_suite_;
		f.new_jump( jump_start, jump_end, cutpoint );

		Size const which_jump = f.jump_nr( jump_start, jump_end );
		f.set_jump_atoms( which_jump,
											jump_start,
											default_jump_atom( pose.residue( jump_start ) ),
											jump_end,
											default_jump_atom( pose.residue( jump_end   ) ) );

		pose.fold_tree( f );

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////
	void
	TransientCutpointHandler::take_out_cutpoints( core::pose::Pose & viewer_pose ){

		using namespace core::chemical;

		Pose pose = viewer_pose; // prevent some conflicts with graphics. Note potential slowdown.

		// remove chainbreak variants. along with fold_tree restorer, put into separate function.
		remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, cutpoint_suite_   );
		remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, cutpoint_suite_ + 1 );

		// return to simple fold tree
		pose.fold_tree( fold_tree_save_ );

		viewer_pose = pose;
	}

} //mover
} //monte_carlo
} //stepwise
} //protocols
