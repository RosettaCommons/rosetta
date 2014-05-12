// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/ResampleMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>
#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/sampling/StepWiseModeler.hh>
#include <protocols/stepwise/sampling/packer/util.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>


using namespace protocols::stepwise::sampling;

static numeric::random::RandomGenerator RG(239145021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.stepwise.monte_carlo.ResampleMover" );


//////////////////////////////////////////////////////////////////////////
// Makes a choice, based on current pose, and information in full_model_info
//  as to where to resample nucleotides or chunks.
//
// Does not add or delete  [see AddOrDeleteMover for those functionalities].
//
//////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	//Constructor
	ResampleMover::ResampleMover(	protocols::stepwise::sampling::StepWiseModelerOP stepwise_modeler ):
		stepwise_modeler_( stepwise_modeler ),
		swa_move_selector_( new SWA_MoveSelector ),
		options_( new StepWiseMonteCarloOptions ),
		minimize_single_res_( false ),
		slide_intermolecular_jumps_( true )
	{}

	//Destructor
	ResampleMover::~ResampleMover()
	{}

	////////////////////////////////////////////////////////////////////
	std::string
	ResampleMover::get_name() const {
		return "ResampleMover";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	ResampleMover::apply( pose::Pose & pose ){
		std::string move_type;
		apply( pose, move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// This version randomly chooses the move element.
	bool
	ResampleMover::apply( pose::Pose & pose,
												std::string & move_type ){

		utility::vector1< SWA_Move > swa_moves;
		swa_move_selector_->set_allow_internal_hinge( options_->allow_internal_hinge_moves() );
		swa_move_selector_->set_allow_internal_local( options_->allow_internal_local_moves() );
		swa_move_selector_->set_intermolecular_frequency( options_->intermolecular_frequency() );
		swa_move_selector_->get_resample_move_elements( pose, swa_moves );

		if ( swa_moves.size() == 0 ) return false;
		SWA_Move const & swa_move = RG.random_element( swa_moves );

		return apply( pose, swa_move, move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ResampleMover::apply( pose::Pose & pose,
												SWA_Move const & swa_move ){
		std::string dummy_move_type;
		return apply( pose, swa_move, dummy_move_type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	ResampleMover::apply( pose::Pose & pose,
												SWA_Move const & swa_move,
												std::string & move_type ){

		using namespace protocols::stepwise;
		using namespace protocols::stepwise::sampling::rna;
		using namespace protocols::stepwise::monte_carlo;
		using namespace core::pose::full_model_info;

		TR << "About to remodel move_element " << swa_move << std::endl;
		move_type = to_string( swa_move.move_type() );
		std::transform(move_type.begin(), move_type.end(), move_type.begin(), ::tolower); // this is why we love C

		// What needs to be set in StepwiseModeler:
		Size remodel_res( 0 ), remodel_suite( 0 ), cutpoint_suite( 0 );

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		Size const num_attachments = swa_move.attachments().size();
		bool const is_single_attachment = ( num_attachments == 1 );
		Size const move_element_size = swa_move.move_element().size();

		if ( is_single_attachment ) {
			remodel_res = get_remodel_res( swa_move, pose );
			if ( slide_intermolecular_jumps_ && swa_move.attachment_type() == JUMP_INTERCHAIN ) slide_jump_randomly( pose, remodel_res );
		} else { // an internal residue or move_element, with two attachments.
			runtime_assert( num_attachments == 2 );
			runtime_assert( swa_move.attachments()[1].attachment_type() == BOND_TO_PREVIOUS );
			runtime_assert( swa_move.attachments()[2].attachment_type() == BOND_TO_NEXT );
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

		bool did_mutation( false );
		// based on 'n' in full_model_info.full_sequence
		if ( move_element_size == 1 ) did_mutation = mutate_res_if_allowed( pose, full_to_sub( swa_move.moving_res(), pose ) );
		bool just_min_after_mutation = ( did_mutation && ( RG.uniform() < options_->just_min_after_mutation_frequency() ) );

		if ( just_min_after_mutation ) {
			stepwise_modeler_->set_moving_res_and_reset( 0 );
			stepwise_modeler_->set_working_prepack_res( packer::figure_out_working_interface_res( pose, remodel_res ) );
		} else {
			stepwise_modeler_->set_moving_res_and_reset( remodel_res );
			stepwise_modeler_->set_figure_out_prepack_res( true );
		}
		utility::vector1< Size > const & moving_res = get_moving_res_from_full_model_info( pose );
		if ( ! minimize_single_res_ ) stepwise_modeler_->set_working_minimize_res( moving_res );

		if ( is_single_attachment ){
			stepwise_modeler_->apply( pose );
		} else {
			runtime_assert( options_->allow_internal_local_moves() );
			TR << "Going to set up TRANSIENT_CUTPOINT_HANDLER with " << remodel_suite << " " << cutpoint_suite << std::endl;
			TransientCutpointHandler cutpoint_handler( remodel_suite, cutpoint_suite );
			if ( ! minimize_single_res_ ) cutpoint_handler.set_minimize_res( moving_res );

			cutpoint_handler.put_in_cutpoints( pose );
			stepwise_modeler_->apply( pose );
			cutpoint_handler.take_out_cutpoints( pose );
		}

		if ( did_mutation ) move_type += "-mut";
		if ( just_min_after_mutation ) move_type = "mut";

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	ResampleMover::get_remodel_res( SWA_Move const & swa_move, pose::Pose const & pose ) const {
		using namespace core::pose::full_model_info;
		runtime_assert( swa_move.attachments().size() == 1 );

		// remodel res will be the first residue in the moving element that is immediatly downstream of the attachment residue.
		Size remodel_res( 0 );
		MoveElement const & move_element = swa_move.move_element();
		for ( Size n = 1; n <= move_element.size(); n++ ){
			Size const & moving_res = full_to_sub( move_element[ n ], pose );
			Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
			if ( parent_res > 0 && sub_to_full( parent_res, pose ) == swa_move.attached_res() ){
				remodel_res = moving_res; break;
			}
		}

		// we may have to reroot pose -- the attachment point might be 'downstream' of the moving element.
		// the rerooting will actually occur later in the Modeler.
		if ( remodel_res == 0 ){
			Size const & moving_res = full_to_sub( swa_move.attached_res(), pose );
			if ( move_element.has_value( sub_to_full( pose.fold_tree().get_parent_residue( moving_res ), pose ) ) ){
				remodel_res = moving_res;
			}
		}

		runtime_assert( remodel_res > 0 );
		return remodel_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	filter_for_proximity( pose::Pose const & pose,
												utility::vector1< Size > & partition_res,
												Size const center_res ) {
		using namespace core::chemical::rna;
		runtime_assert( partition_res.has_value( center_res ) );
		static Distance const proximity_cutoff( 8.0 );
		utility::vector1< Size > filtered_partition_res;
		Vector const & center_xyz = pose.residue( center_res ).xyz( default_jump_atom( pose.residue( center_res ) ) );
		for ( Size n = 1; n <= partition_res.size(); n++ ){
			Size const new_res = partition_res[ n ];
			Vector const & new_xyz =  pose.residue( new_res ).xyz( default_jump_atom( pose.residue( new_res ) ) );
			if ( ( new_xyz - center_xyz ).length() < proximity_cutoff ) filtered_partition_res.push_back( new_res );
		}
		partition_res = filtered_partition_res;
		runtime_assert( partition_res.has_value( center_res ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	ResampleMover::slide_jump_randomly( pose::Pose & pose, Size & remodel_res ) const {

		using namespace core::kinematics;
		using namespace core::chemical::rna;

		FoldTree f = pose.fold_tree();
		Size const jump_nr = f.get_jump_that_builds_residue( remodel_res );
		Size const reference_res = f.upstream_jump_residue( jump_nr );

		utility::vector1< Size > root_partition_res, moving_partition_res;
		figure_out_root_and_moving_partition_res( pose, remodel_res, root_partition_res, moving_partition_res );

		if ( options_->local_redock_only() ){
			filter_for_proximity( pose, root_partition_res, reference_res );
			filter_for_proximity( pose, moving_partition_res, remodel_res );
		}

		// need to make sure JUMP_INTERMOL remain in different chains!
		utility::vector1< Size > chains = 	figure_out_chains_from_full_model_info_const( pose );
		utility::vector1< std::pair< Size, Size > > possible_jump_pairs;
		for ( Size i = 1; i <= root_partition_res.size(); i++ ) {
			Size const & root_res = root_partition_res[ i ];
			for ( Size j = 1; j <= moving_partition_res.size(); j++ ) {
				Size const & move_res = moving_partition_res[ j ];
				if ( chains[ root_res ] != chains[ move_res ] ) possible_jump_pairs.push_back( std::make_pair( root_res, move_res ) );
			}
		}
		std::pair< Size, Size > const new_jump_pair = RG.random_element( possible_jump_pairs );
		Size const new_reference_res = new_jump_pair.first;
		Size const new_remodel_res   = new_jump_pair.second;

		f.slide_jump( jump_nr, new_reference_res, new_remodel_res );
		f.set_jump_atoms( jump_nr, default_jump_atom( pose.residue( new_reference_res ) ),
											default_jump_atom( pose.residue( new_remodel_res ) ) );

		pose.fold_tree( f );
		TR << "Slid jump from: " << remodel_res << "--" << reference_res <<
			" to: " << new_remodel_res << " -- " << new_reference_res << std::endl;

		remodel_res = new_remodel_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	ResampleMover::set_options( StepWiseMonteCarloOptionsCOP options ){
		options_ = options;
	}


} //rna
} //monte_carlo
} //stepwise
} //protocols
