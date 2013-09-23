// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/RNA_ResampleMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/monte_carlo/RNA_ResampleMover.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/TransientCutpointHandler.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>
#include <protocols/swa/monte_carlo/SWA_MonteCarloUtil.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>


static numeric::random::RandomGenerator RG(239145021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.swa.monte_carlo.RNA_ResampleMover" );

namespace protocols {
namespace swa {
namespace monte_carlo {

	//Constructor
	RNA_ResampleMover::RNA_ResampleMover(	protocols::swa::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler ):
		stepwise_rna_modeler_( stepwise_rna_modeler ),
		just_min_after_mutation_frequency_( 0.5 ),
		allow_internal_moves_( false ),
		minimize_single_res_( false )
	{}

	//Destructor
	RNA_ResampleMover::~RNA_ResampleMover()
	{}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This version randomly chooses the move element.
	bool
	RNA_ResampleMover::apply( pose::Pose & pose,
														std::string & move_type ){

		utility::vector1< SWA_Move > swa_moves;
		get_resample_terminal_move_elements( pose, swa_moves );
		if ( allow_internal_moves_ ) get_resample_internal_move_elements( pose, swa_moves );

		//		TR << "POSSIBLE RESAMPLE MOVES! ";
		//		for ( Size n = 1; n <= swa_moves.size(); n++ ) TR << swa_moves[n];
		//		TR << std::endl;

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

		using namespace protocols::swa;
		using namespace protocols::swa::rna;
		using namespace protocols::swa::monte_carlo;
		using namespace core::pose::full_model_info;

		TR << "About to remodel move_element " << swa_move << std::endl;

		// What needs to be set in StepwiseRNA_Modeler:
		Size remodel_res( 0 ), remodel_suite( 0 ), cutpoint_suite( 0 );

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		Size const num_attachments = swa_move.attachments().size();
		bool const is_at_terminus = ( num_attachments == 1 );
		Size const move_element_size = swa_move.move_element().size();

		if ( is_at_terminus ) {
			Attachment const & attachment = swa_move.attachments()[ 1 ];
			AttachmentType const & attachment_type = attachment.attachment_type();
			if ( move_element_size == 1) {
				remodel_res = res_list.index( swa_move.moving_res() );
				remodel_suite = ( attachment_type == ATTACHED_TO_PREVIOUS ) ? remodel_res - 1 : remodel_res;
			} else { // remodel res will be internal... need to supply *suite* number.
				Size const & attachment_res = attachment.attached_res();
				if ( attachment_type == ATTACHED_TO_PREVIOUS ){
					remodel_res = res_list.index( attachment_res );
				} else {
					runtime_assert( attachment_type == ATTACHED_TO_NEXT );
					remodel_res = res_list.index( attachment_res ) - 1;
					runtime_assert( remodel_res == res_list.index( attachment_res - 1 ) );
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
				runtime_assert( cutpoint_suite == res_list.index( attached_res_next - 1 ) );
			}
		}
		runtime_assert( remodel_res > 0 );

		bool did_mutation( false );
		// based on 'n' in full_model_info.full_sequence
		if ( move_element_size == 1 ) did_mutation = mutate_res_if_allowed( pose, full_to_sub( swa_move.moving_res(), pose ) );
		bool just_min_after_mutation_ = ( did_mutation && ( RG.uniform() < just_min_after_mutation_frequency_ ) );

		stepwise_rna_modeler_->set_moving_res_and_reset( remodel_res );
		stepwise_rna_modeler_->set_skip_sampling( just_min_after_mutation_ );

		// LATER SHOULD REPLACE THIS WITH FIXED DOMAIN MAP -- NEED TO UPDATE STEPWISE MODELER -- rhiju.
		utility::vector1< Size > const & moving_res = get_moving_res_from_full_model_info( pose );
		if ( ! minimize_single_res_ ) stepwise_rna_modeler_->set_minimize_res( moving_res );

		if ( is_at_terminus ){

			move_type = "resample_terminus";
			stepwise_rna_modeler_->apply( pose );

		} else {
			runtime_assert( allow_internal_moves_ );
			move_type = "resample_internal_local";

			TransientCutpointHandler cutpoint_handler( remodel_suite, cutpoint_suite );
			if ( ! minimize_single_res_ ) cutpoint_handler.set_minimize_res( moving_res );

			cutpoint_handler.put_in_cutpoints( pose );
			stepwise_rna_modeler_->apply( pose );
			cutpoint_handler.take_out_cutpoints( pose );
		}

		if ( did_mutation ) move_type += "-mut";
		if ( just_min_after_mutation_ ) move_type = "mut";

		return true;

	}


} //monte_carlo
} //swa
} //protocols
