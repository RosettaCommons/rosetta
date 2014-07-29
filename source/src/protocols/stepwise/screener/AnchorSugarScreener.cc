// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/AnchorSugarScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/AnchorSugarScreener.hh>
#include <protocols/stepwise/screener/AnchorSugarScreener.fwd.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/moves/CompositionMover.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.AnchorSugarScreener" );

using namespace protocols::stepwise::modeler::rna::sugar;
using namespace protocols::stepwise::modeler::rna::checker;
using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	AnchorSugarScreener::AnchorSugarScreener( SugarModeling const & anchor_sugar_modeling,
																						RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_to_anchor_checker,
																						pose::Pose & sugar_screening_pose,
																						bool const is_prepend,
																						RNA_AtrRepCheckerOP atr_rep_checker_with_instantiated_sugar,
																						utility::vector1< RNA_AtrRepCheckerOP > atr_rep_checkers_for_anchor_sugar_models,
																						TagDefinitionOP tag_definition ):
		anchor_sugar_modeling_( anchor_sugar_modeling ),
		chain_closable_geometry_to_anchor_checker_( chain_closable_geometry_to_anchor_checker ),
		sugar_screening_pose_( sugar_screening_pose ),
		atr_rep_checker_with_instantiated_sugar_( atr_rep_checker_with_instantiated_sugar ),
		atr_rep_checkers_for_anchor_sugar_models_( atr_rep_checkers_for_anchor_sugar_models ),
		tag_definition_( tag_definition ),
		is_prepend_( is_prepend ),
		moving_atom_name_( ( is_prepend ) ? " O3'" : " C5'" ),
		reference_atom_name_( ( is_prepend ) ? " C5'" : " O3'" ),
		anchor_sugar_solution_number_( 0 )
	{}

	//Destructor
	AnchorSugarScreener::~AnchorSugarScreener()
	{}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// following function might be dramatically simplifiable, in which case we may want to return
	// it to main loop.
	bool
	AnchorSugarScreener::check_screen(){

		anchor_sugar_solution_number_ = 0;

		//OK check that with this sugar, the chain can be theoretically closed..
		if ( !anchor_sugar_modeling_.sample_sugar ){

			bool const ok = chain_closable_geometry_to_anchor_checker_->check_screen( sugar_screening_pose_ );
			//			TR << "DIST_SQUARED " << chain_closable_geometry_to_anchor_checker_->dist_squared() << " " << ok << std::endl;
			if ( !ok ) return 0;

			if ( atr_rep_checker_with_instantiated_sugar_ &&
					 !atr_rep_checker_with_instantiated_sugar_->check_screen( sugar_screening_pose_ ) ) return false; // wait a minute... why is this in here? oh, because base can move in different sampled sugar modeling conformations.

			return true;
		}

		// The point of this section is to look for *any* conformation of sugar in anchor residue that passes
		// screens, going from the lowest energy option on up.
		//Ok, since anchor_sugar_modeling_.pose_list is sorted by SCORE, the lower energy conformations are tried first!
		for ( Size n = 1; n <= anchor_sugar_modeling_.pose_list.size(); n++ ){
			pose::Pose const & anchor_sugar_modeling_pose = *anchor_sugar_modeling_.pose_list[n];

			// THIS WAS WORKING -- KEEP THIS IN
			// is_prepend --> sugar_screening_pose [moving] = 5' pose [need O3'], anchor_sugar = 3' pose [need C5']
			if ( !chain_closable_geometry_to_anchor_checker_->check_screen( sugar_screening_pose_, anchor_sugar_modeling_pose, is_prepend_) ) continue;

			// following could be replaced with (pre-instantiated) CopyDofMover -- see below.
			copy_bulge_res_and_sugar_torsion( anchor_sugar_modeling_, sugar_screening_pose_, anchor_sugar_modeling_pose );


			// DO NOT CHECK IN
			//			if ( !chain_closable_geometry_to_anchor_checker_->check_screen( sugar_screening_pose_ ) ) continue;

			// DO NOT CHECK IN.
			//			if ( atr_rep_checker_with_instantiated_sugar_ &&
			//					 !atr_rep_checker_with_instantiated_sugar_->check_screen( sugar_screening_pose_ ) ) continue;

			// THIS IS THE RIGHT THING TO DO.
			// This is in here because the anchor sugar models can have slightly shifted bases (not just riboses!) due to a minimization step that
			// can occur in VirtualRiboseSampler [see the option: do_minimize].
			if ( atr_rep_checkers_for_anchor_sugar_models_[ n ] &&
					 !atr_rep_checkers_for_anchor_sugar_models_[ n ]->check_screen( sugar_screening_pose_ ) ) continue;

			anchor_sugar_solution_number_ = n;

			tag_definition_->append_to_tag( tag_from_pose( anchor_sugar_modeling_pose ) );
			return true;
		}

		return false;
	}

/////////////////////////////////////////////////////////
void
AnchorSugarScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
	if ( !anchor_sugar_modeling_.sample_sugar ){
		update_mover->add_mover( 0 );
		restore_mover->add_mover( 0 );
		return;
	}

	runtime_assert( anchor_sugar_solution_number_ > 0  );
	Size const n = anchor_sugar_solution_number_;
	pose::Pose const & anchor_sugar_modeling_pose = *anchor_sugar_modeling_.pose_list[n];

	std::map< core::Size, core::Size > res_map;  //This is map from sub numbering to input_res numbering..
	res_map[ anchor_sugar_modeling_.moving_res    ] = anchor_sugar_modeling_.moving_res;
	if ( anchor_sugar_modeling_.bulge_res > 0 ) res_map[ anchor_sugar_modeling_.bulge_res     ] = anchor_sugar_modeling_.bulge_res;
	res_map[ anchor_sugar_modeling_.reference_res ] = anchor_sugar_modeling_.reference_res;

	simple_moves::CopyDofMoverOP copy_dof_mover = new simple_moves::CopyDofMover( anchor_sugar_modeling_pose, res_map );
	update_mover->add_mover( copy_dof_mover );
	restore_mover->add_mover( 0 );
}


} //screener
} //stepwise
} //protocols
