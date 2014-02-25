// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/rna/RNA_ResampleMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/rna/RNA_ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>
#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.hh>
#include <protocols/stepwise/monte_carlo/rna/TransientCutpointHandler.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>


static numeric::random::RandomGenerator RG(239145021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.stepwise.monte_carlo.RNA_ResampleMover" );


//////////////////////////////////////////////////////////////////////////
// Makes a choice, based on current pose, and information in full_model_info
//  as to where to resample nucleotides or chunks.
//
// Does not add or delete  [see RNA_AddOrDeleteMover for those functionalities].
//
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	//Constructor
	RNA_ResampleMover::RNA_ResampleMover(	protocols::stepwise::sampling::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler ):
		stepwise_rna_modeler_( stepwise_rna_modeler ),
		swa_move_selector_( new SWA_MoveSelector ),
		options_( new StepWiseRNA_MonteCarloOptions ),
		minimize_single_res_( false )
	{}

	//Destructor
	RNA_ResampleMover::~RNA_ResampleMover()
	{}

	////////////////////////////////////////////////////////////////////
	std::string
	RNA_ResampleMover::get_name() const {
		return "RNA_ResampleMover";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ResampleMover::apply( pose::Pose & pose ){
		std::string move_type;
		apply( pose, move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This version randomly chooses the move element.
	bool
	RNA_ResampleMover::apply( pose::Pose & pose,
														std::string & move_type ){

		utility::vector1< SWA_Move > swa_moves;
		swa_move_selector_->set_allow_internal_hinge( options_->allow_internal_hinge_moves() );
		swa_move_selector_->set_allow_internal_local( options_->allow_internal_local_moves() );
		swa_move_selector_->get_resample_move_elements( pose, swa_moves );

		TR.Debug << "POSSIBLE RESAMPLE MOVES!";
		for ( Size n = 1; n <= swa_moves.size(); n++ ) TR.Debug << " " << swa_moves[n];
		TR.Debug << std::endl;

		if ( swa_moves.size() == 0 ) return false;
		SWA_Move const & swa_move = RG.random_element( swa_moves );

		return apply( pose, swa_move, move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ResampleMover::apply( pose::Pose & pose,
														SWA_Move & swa_move ){
		std::string dummy_move_type;
		return apply( pose, swa_move, dummy_move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ResampleMover::apply( pose::Pose & pose,
														SWA_Move const & swa_move,
														std::string & move_type ){

		using namespace protocols::stepwise;
		using namespace protocols::stepwise::sampling::rna;
		using namespace protocols::stepwise::monte_carlo;
		using namespace core::pose::full_model_info;

		TR << "About to remodel move_element " << swa_move << std::endl;

		// What needs to be set in StepwiseRNA_Modeler:
		Size remodel_res( 0 ), remodel_suite( 0 ), cutpoint_suite( 0 );

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		Size const num_attachments = swa_move.attachments().size();
		bool const is_single_attachment = ( num_attachments == 1 );
		Size const move_element_size = swa_move.move_element().size();

		if ( is_single_attachment ) {
			Attachment const & attachment = swa_move.attachments()[ 1 ];
			AttachmentType const & attachment_type = attachment.attachment_type();
			//			TR << TR.Red << "Attachment " << attachment << TR.Reset << std::endl;
			if ( move_element_size == 1) {
				remodel_res = res_list.index( swa_move.moving_res() );
				remodel_suite = ( attachment_type == ATTACHED_TO_PREVIOUS ) ? remodel_res - 1 : remodel_res;
				// also possible that this is a 'floating base'
				if ( attachment_type == JUMP_TO_PREV_IN_CHAIN || attachment_type == JUMP_TO_NEXT_IN_CHAIN ) remodel_suite = 0;
			} else { // remodel res will be internal... need to supply *suite* number.
				Size const & attachment_res = attachment.attached_res();
				if ( attachment_type == ATTACHED_TO_PREVIOUS ){
					remodel_res = res_list.index( attachment_res );
				} else {
					runtime_assert( attachment_type == ATTACHED_TO_NEXT ); // cannot yet handle jumps.
					remodel_res = res_list.index( attachment_res ) - 1;
					runtime_assert( static_cast<int>(remodel_res) == res_list.index( attachment_res - 1 ) );
				}
				remodel_suite = remodel_res;
			}
		} else { // an internal residue or move_element, with two attachments.
			runtime_assert( num_attachments == 2 );
			runtime_assert( swa_move.attachments()[1].attachment_type() == ATTACHED_TO_PREVIOUS );
			runtime_assert( swa_move.attachments()[2].attachment_type() == ATTACHED_TO_NEXT );
			if ( move_element_size == 1) { // single residue
				remodel_res   = res_list.index( swa_move.moving_res() );
				remodel_suite = remodel_res - 1;
				cutpoint_suite  = res_list.index( swa_move.moving_res() );
			} else { // fixed move_element
				Size const attached_res_prev = swa_move.attachments()[1].attached_res();
				Size const attached_res_next = swa_move.attachments()[2].attached_res();
				remodel_res   = res_list.index( attached_res_prev ); // this is now the *suite*
				remodel_suite = remodel_res;
				cutpoint_suite  = res_list.index( attached_res_next ) - 1;
				runtime_assert( static_cast<int>(cutpoint_suite) == res_list.index( attached_res_next - 1 ) );
			}
		}
		runtime_assert( remodel_res > 0 );

		// trying to catch last bugs.
		TR << pose.fold_tree() << TR.Reset;
		TR << pose.annotated_sequence() << std::endl;
		//		TR << TR.Red << "Is at terminus: " << is_single_attachment << " num attachments " << num_attachments << TR.Reset << std::endl;

		bool did_mutation( false );
		// based on 'n' in full_model_info.full_sequence
		if ( move_element_size == 1 ) did_mutation = mutate_res_if_allowed( pose, full_to_sub( swa_move.moving_res(), pose ) );
		bool just_min_after_mutation = ( did_mutation && ( RG.uniform() < options_->just_min_after_mutation_frequency() ) );

		stepwise_rna_modeler_->set_moving_res_and_reset( remodel_res );
		stepwise_rna_modeler_->set_skip_sampling( just_min_after_mutation );

		// LATER SHOULD REPLACE THIS WITH FIXED DOMAIN MAP -- NEED TO UPDATE STEPWISE MODELER -- rhiju.
		utility::vector1< Size > const & moving_res = get_moving_res_from_full_model_info( pose );
		if ( ! minimize_single_res_ ) stepwise_rna_modeler_->set_minimize_res( moving_res );

		clear_constraints_recursively( pose );
		if ( get_native_pose() ) {
			superimpose_recursively_and_add_constraints( pose, *get_native_pose(),
																									 options_->constraint_x0(), options_->constraint_tol() );
		}

		if ( is_single_attachment ){

			move_type = "resample_terminus";
			stepwise_rna_modeler_->apply( pose );

		} else {
			runtime_assert( options_->allow_internal_local_moves() );
			move_type = "resample_internal_local";

			TR << "Going to set up TRANSIENT_CUTPOINT_HANDLER with " << remodel_suite << " " << cutpoint_suite << std::endl;
			TransientCutpointHandler cutpoint_handler( remodel_suite, cutpoint_suite );
			if ( ! minimize_single_res_ ) cutpoint_handler.set_minimize_res( moving_res );

			cutpoint_handler.put_in_cutpoints( pose );
			stepwise_rna_modeler_->apply( pose );
			cutpoint_handler.take_out_cutpoints( pose );
		}

		if ( did_mutation ) move_type += "-mut";
		if ( just_min_after_mutation ) move_type = "mut";

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ResampleMover::set_options( StepWiseRNA_MonteCarloOptionsCOP options ){
		options_ = options;
	}


} //rna
} //monte_carlo
} //stepwise
} //protocols
