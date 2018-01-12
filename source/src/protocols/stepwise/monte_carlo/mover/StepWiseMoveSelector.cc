// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/mover/options/StepWiseMoveSelectorOptions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh> // probably do not need RNA in here.
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/util.hh> // Probably could get rid of RNA-specialized stuff here.
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/SubMotifInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/loop_graph/LoopGraph.hh>
#include <utility>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>
#include <utility/vector1.functions.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

#include <map>
#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.stepwise.monte_carlo.StepWiseMoveSelector" );

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;
using namespace protocols::stepwise::modeler;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Now with detailed balance! See how this is called in stepwise/moves/StepWiseMasterMover,
//   including use of proposal_probabilities.
//
// Historically, this class accreted out of some simple functions to figure out what could
//   be added or deleted to a pose; probably could be rewritten now to make logic
//   clear. In particular, 'delete' and 'resample' both involve basicaly the same
//   subset of internal DOFs already existing in a pose -- could be computed at the
//   same time. Also, less emphasis now on moves for a specific pose, as opposed to
//   all poses involved in a full model (encoded in other_pose_list in FullModelInfo) -- this works in
//   new_move_selector mode (becoming default in late 2015).
//
// dock_domain defines two domains that might be docked. If they are alread docked, don't do any moves
//  that might dock them 'again'. This data structure is rather limiting to two dock domains, and was originally
//  introduced to avoid 'complex loop cycles' (nested cycles), but the code can now filter for those separately.
//  we might want to deprecate this check.
//
//  Also, StepWiseMove should probably be refactored -- it really just requires 2 residue sets that
//   define partitions connected by the DOF to add, delete, or resample. The order
//   of the partitions is no longer quite relevant (and can change during modeling).
//
//   But in ADD moves, we just specify one of the added residues.
//   And in all moves, we just specify one of the 'attached' residues. We do need to specify attachment
//    positions when it is ambiguous (e.g., filling in a single nt that could be attached on either side.)
//
// New StepWiseMove should have:
//    partition1, partition2 (order does not matter, but use convention where partition1 is the one with lower residue number.)
//    connections [not ordered: res in partition 1, res in partition2, BOND/JUMP]
//
// Several of the core routines could be accelerated, if needed. E.g., domain_map
//   can be computed once; and instead of using sample_res.has_value(...), we
//   could access lookup maps.
//
//  -- rhiju, jan 2015.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

//Constructor
StepWiseMoveSelector::StepWiseMoveSelector( options::StepWiseMoveSelectorOptionsCOP options ):
	options_(std::move( options )),
	allow_delete_( true ),
	force_unique_moves_( false ),
	choose_random_( true )
{}

//Constructor
StepWiseMoveSelector::StepWiseMoveSelector():
	options_( options::StepWiseMoveSelectorOptionsCOP( new options::StepWiseMoveSelectorOptions ) ),
	allow_delete_( true ),
	force_unique_moves_( false ),
	choose_random_( true )
{
}

//Destructor
StepWiseMoveSelector::~StepWiseMoveSelector() = default;

/// @brief copy constructor
StepWiseMoveSelector::StepWiseMoveSelector( StepWiseMoveSelector const & src ) :
	ReferenceCount( src )
{
	*this = src;
}

/// @brief clone the options
StepWiseMoveSelectorOP
StepWiseMoveSelector::clone() const
{
	return StepWiseMoveSelectorOP( new StepWiseMoveSelector( *this ) );
}


bool
StepWiseMoveSelector::figure_out_all_possible_moves( pose::Pose const & pose,
	bool const verbose /* = true */ ) {

	swa_moves_.clear();
	proposal_probabilities_.clear();

	// resampling, folding, and docking.
	// find moves for pose on which we have 'focus'.
	fill_moves_for_pose( pose );

	// find moves for other poses -- but weight those less
	// this must occur after fill_moves_for_pose, and before denovo_moves
	fill_moves_for_other_poses( pose );

	// de novo stuff. independent of the specific pose that we're looking at.
	fill_denovo_moves( pose );

	// in design runs, can vary lengths of designable loops, in or connecting the different poses
	fill_vary_loop_length_moves( pose );

	normalize_probabilities();
	if ( verbose ) output_moves();

	return ( swa_moves_.size() > 0 );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::output_moves() const {
	runtime_assert( swa_moves_.size() == proposal_probabilities_.size() );
	for ( Size n = 1; n <= swa_moves_.size(); n++ ) {
		TR << TR.Green << ObjexxFCL::format::F( 8, 3, proposal_probabilities_[n] ) << ": " << swa_moves_[ n ] << TR.Reset << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::fill_moves_for_pose( pose::Pose const & pose ) {

	utility::vector1< StepWiseMove > resample_moves, intramolecular_add_and_delete_moves, docking_add_and_delete_moves;

	// combining "add_and_delete" into one function was historical --
	// might be better to have delete moves paired with resample moves, as they both look at
	// existing, *internal*  DOFs. And have add moves separated, as they involve new DOFs.
	get_intramolecular_add_and_delete_elements( pose, intramolecular_add_and_delete_moves );

	save_moves( intramolecular_add_and_delete_moves, options_->add_delete_frequency()  );

	get_docking_add_and_delete_elements( pose, docking_add_and_delete_moves );
	save_moves( docking_add_and_delete_moves, options_->docking_frequency()  );

	get_resample_move_elements( pose, resample_moves, false /* save_moves */ );
	save_moves( resample_moves, ( 1 - options_->add_delete_frequency() )  );
}


/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::fill_denovo_moves( pose::Pose const & pose ) {
	utility::vector1< StepWiseMove >  from_scratch_add_and_delete_moves, submotif_add_moves;

	get_from_scratch_add_and_delete_elements( pose, from_scratch_add_and_delete_moves );
	save_moves( from_scratch_add_and_delete_moves, options_->from_scratch_frequency()  );

	get_submotif_add_moves( pose, submotif_add_moves );
	save_moves( submotif_add_moves, options_->submotif_frequency() );

	// not implemented yet ...
	// get_submotif_seed_add_moves( pose, submotif_add_moves );
	// save_moves( submotif_seed_add_moves, submotif_seed_frequency() );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::fill_vary_loop_length_moves( pose::Pose const & pose ) {
	if ( options_->vary_loop_length_frequency() == 0 ) return;

	utility::vector1< StepWiseMove >  vary_loop_length_moves;
	get_vary_loop_length_moves( pose, vary_loop_length_moves );

	save_moves( vary_loop_length_moves, options_->vary_loop_length_frequency() );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// make sure that we get back to original pose...
// and set weights of other moves so that their sum is some specific fraction of
// the original weight.
void
StepWiseMoveSelector::fill_moves_for_other_poses( pose::Pose const & pose ) {

	Size const move_idx_original = swa_moves_.size();
	Real const original_weight = sum_probabilities( 1, move_idx_original );

	// use a residue in each pose to serve as a unique identifier of each pose.
	utility::vector1< Size > first_res_in_other_poses;
	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ) {
		first_res_in_other_poses.push_back( const_full_model_info( *other_pose_list[ n ] ).res_list()[ 1 ] );
	}

	Pose pose_to_switch = pose;
	for ( Size n = 1; n <= first_res_in_other_poses.size(); n++ ) {
		Size const other_pose_idx = const_full_model_info( pose_to_switch ).get_idx_for_other_pose_with_residue( first_res_in_other_poses[ n ] );
		switch_focus_to_other_pose( pose_to_switch, other_pose_idx );
		fill_moves_for_pose( pose_to_switch );
	}

	Size const move_idx_final = swa_moves_.size();
	normalize_probabilities( move_idx_original + 1, move_idx_final,
		original_weight * options_->switch_focus_frequency() );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// save these moves so that the sum of the weights comes out to 'total_weight'.
void
StepWiseMoveSelector::save_moves( utility::vector1< StepWiseMove > const & moves,
	Real const total_weight,
	bool const clear_moves_before_saving /* = false */) {
	if ( clear_moves_before_saving ) {
		swa_moves_.clear();
		proposal_probabilities_.clear();
	}

	if ( total_weight == 0.0 ) return;

	utility::vector1< Real > relative_weights;
	for ( Size n = 1; n <= moves.size(); n++ ) {
		// a bit hacky here -- reduce relative weight of skip_bulge moves.
		if ( ( moves[ n ].move_type() == ADD ||
				moves[ n ].move_type() == ADD_SUBMOTIF ) &&
				( moves[ n ].attachment_type() == JUMP_TO_PREV_IN_CHAIN ||
				moves[ n ].attachment_type() == JUMP_TO_NEXT_IN_CHAIN ) ) {
			relative_weights.push_back( options_->skip_bulge_frequency() );
		} else {
			relative_weights.push_back( 1 );
		}
	}

	Real relative_weight_sum( 0.0 );
	for ( Size n = 1; n <= moves.size(); n++ ) relative_weight_sum += relative_weights[ n ];

	for ( Size n = 1; n <= moves.size(); n++ ) {
		if ( force_unique_moves_ ) {
			if ( swa_moves_.has_value( moves[ n ] ) ) TR << moves[ n ] <<  " is already in " << swa_moves_ << std::endl;
			runtime_assert( !swa_moves_.has_value( moves[ n ] ) );
		}
		swa_moves_.push_back( moves[ n ] );
		proposal_probabilities_.push_back( total_weight * relative_weights[ n ]/ relative_weight_sum );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseMoveSelector::sum_probabilities( Size const start_idx /* = 1 */,
	Size const final_idx_input /* = 0 */ ) const
{
	Size const final_idx = ( final_idx_input > 0 ) ? final_idx_input : proposal_probabilities_.size();
	Real weight( 0.0 );
	for ( Size n = start_idx; n <= final_idx; n++ ) {
		weight += proposal_probabilities_[ n ];
	}
	return weight;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::normalize_probabilities( Size const start_idx /* = 1 */,
	Size const final_idx_input /* = 0 */,
	Real const desired_weight /* = 1 */ ) {
	Size const final_idx = ( final_idx_input > 0 ) ? final_idx_input : proposal_probabilities_.size();
	Real weight = sum_probabilities( start_idx, final_idx );
	for ( Size n = start_idx; n <= final_idx; n++ ) {
		proposal_probabilities_[ n ] *= desired_weight / weight;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// return a move, weighted in probability by proposal_probabilities.
StepWiseMove
StepWiseMoveSelector::select_random_move( pose::Pose const &  ) const {
	runtime_assert( swa_moves_.size() > 0 );
	runtime_assert( swa_moves_.size() == proposal_probabilities_.size() );
	Real const r = numeric::random::rg().uniform();
	Real weight( 0.0 );
	for ( Size n = 1; n <= proposal_probabilities_.size(); n++ ) {
		weight += proposal_probabilities_[ n ];
		if ( r < weight ) return swa_moves_[ n ];
	}
	runtime_assert( std::abs( weight - 1.0 ) < 1.0e-5 );
	return swa_moves_[ swa_moves_.size() ];
}

/////////////////////////////////////////////////////////////////////////////////////////////////
Real
StepWiseMoveSelector::proposal_probability( StepWiseMove const & swa_move ) const {
	runtime_assert( swa_moves_.size() == proposal_probabilities_.size() );
	if ( !swa_moves_.has_value( swa_move ) ) {
		TR << TR.Red << "WARNING! WARNING! WARNING! swa_moves_.has_value( swa_move ) == false!" << TR.Reset << std::endl;
		return 1.0;
	}
	runtime_assert( swa_moves_.has_value( swa_move ) );
	return proposal_probabilities_[ swa_moves_.index( swa_move ) ];
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_add_or_delete_element( pose::Pose const & pose,
	StepWiseMove & swa_move ) {

	using namespace core::pose::full_model_info;

	utility::vector1< StepWiseMove > swa_moves;
	if ( choose_random_ ) {
		if ( numeric::random::rg().uniform() < options_->docking_frequency() ) {
			TR << "made it" << std::endl;
			get_docking_add_and_delete_elements( pose, swa_moves ); // docking
		} else if ( numeric::random::rg().uniform() < options_->from_scratch_frequency() ) {
			get_from_scratch_add_and_delete_elements( pose, swa_moves );
		} else if ( numeric::random::rg().uniform() < options_->submotif_frequency() ) {
			get_submotif_add_moves( pose, swa_moves );
			// } else if ( numeric::random::rg().uniform() < options_->vary_loop_length_frequency() ) {
			//  get_vary_loop_length_moves( pose, swa_moves );
		}
	}

	if ( swa_moves.size() == 0  ) get_intramolecular_add_and_delete_elements( pose, swa_moves ); // folding

	// backup -- in case we're starting really from nil.
	if ( swa_moves.size() == 0  ) get_docking_add_and_delete_elements( pose, swa_moves ); // docking
	if ( swa_moves.size() == 0  ) get_from_scratch_add_and_delete_elements( pose, swa_moves ); // folding
	if ( swa_moves.size() == 0  ) get_submotif_add_moves( pose, swa_moves ); // special submotifs
	//if ( swa_moves.size() == 0  ) get_vary_loop_length_moves( pose, swa_moves ); // vary loop length

	save_moves( swa_moves, 1, true /* clear_moves_before_saving */ );

	if ( swa_moves.size() == 0 ) {
		swa_move = StepWiseMove();
		return;
	}

	if ( choose_random_ ) {
		swa_move = numeric::random::rg().random_element( swa_moves );
	} else  {
		swa_move = swa_moves[ 1 ];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_resample_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves,
	bool save_resample_moves /* = true */ ) {
	if ( options_->allow_internal_hinge_moves() ) {
		get_resample_internal_move_elements( pose, swa_moves );
	} else {
		get_resample_terminal_move_elements( pose, swa_moves );
	}

	if ( options_->allow_internal_local_moves() ) {
		get_resample_internal_local_move_elements( pose, swa_moves );
	}

	if ( save_resample_moves ) { // when called as stand-alone.
		swa_moves_.clear();
		save_moves( swa_moves, 1, true /*clear_moves_before_saving*/ );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_resample_terminal_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	utility::vector1< StepWiseMove > swa_moves_terminal;
	get_terminal_move_elements( pose, swa_moves_terminal, RESAMPLE );
	remove_from_consideration_first_multi_residue_move_element( swa_moves_terminal, true /*remove_even_if_not_singlet*/ );
	swa_moves.reserve( swa_moves.size() + swa_moves_terminal.size() );
	std::move( std::begin( swa_moves_terminal ), std::end( swa_moves_terminal ), std::back_inserter( swa_moves ) );
	//for ( auto const & Size n = 1; n <= swa_moves_terminal.size(); n++ ) {
	// swa_moves.push_back( swa_moves_terminal[ n ] );
	//}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_resample_internal_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {
	utility::vector1< StepWiseMove > swa_moves_internal;
	get_intramolecular_split_move_elements( pose, swa_moves_internal, RESAMPLE );
	if ( options_->docking_frequency() > 0.0 ) get_docking_split_move_elements( pose, swa_moves_internal, RESAMPLE );
	swa_moves.reserve( swa_moves.size() + swa_moves_internal.size() );
	std::move( std::begin( swa_moves_internal), std::end( swa_moves_internal ), std::back_inserter( swa_moves ) );
	//for ( Size n = 1; n <= swa_moves_internal.size(); n++ ) {
	// swa_moves.push_back( swa_moves_internal[ n ] );
	//}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_resample_internal_local_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	utility::vector1< StepWiseMove > swa_moves_internal;
	get_internal_move_elements( pose, swa_moves_internal, RESAMPLE_INTERNAL_LOCAL );
	// don't delete a multi_residue_move_element if its the only one!
	remove_from_consideration_first_multi_residue_move_element( swa_moves_internal, false /*remove_even_if_not_singlet*/ );
	swa_moves.reserve( swa_moves.size() + swa_moves_internal.size() );
	std::move( std::begin( swa_moves_internal ), std::end( swa_moves_internal ), std::back_inserter( swa_moves ) );
	//for ( Size n = 1; n <= swa_moves_internal.size(); n++ ) swa_moves.push_back( swa_moves_internal[ n ] );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// intramolecular [in-chain]
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_intramolecular_add_and_delete_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	if ( allow_delete_ )     get_intramolecular_delete_move_elements( pose, swa_moves );
	get_intramolecular_add_move_elements( pose, swa_moves );

	if ( force_unique_moves_ ) filter_pose_order( swa_moves, pose );
	if ( options_->filter_complex_cycles() ) filter_complex_cycles( swa_moves, pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_intramolecular_add_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ){

	using namespace core::pose::full_model_info;

	Size const & nres( pose.size() );
	kinematics::FoldTree const & fold_tree( pose.fold_tree() );
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	Size nres_full = full_model_info.full_sequence().size();
	utility::vector1< StepWiseMove > swa_moves_skip_bulge; // collect separately and then add to swa_moves at end.

	if ( res_list.size() == 0 ) return;

	utility::vector1< bool > is_cutpoint_in_full_pose( nres_full, false );
	for ( Size n = 1; n <= cutpoint_open_in_full_model.size(); n++ ) is_cutpoint_in_full_pose[ cutpoint_open_in_full_model[n] ] = true;

	for ( Size i = 1; i <= nres; i++ ) {

		if ( ( i != nres ) &&
				( !fold_tree.is_cutpoint( i ) || (res_list[ i ]+1 >= res_list[ i+1 ]) ) ) continue;

		// could be a 3' chain terminus
		Size const i_full = res_list[ i ] ;
		if ( i_full >= nres_full || is_cutpoint_in_full_pose[ i_full ] ) continue;

		// direct addition of single residue (or potentially, a block of residues starting with this residue)
		if ( is_addable_res( i_full + 1, pose ) ) {
			swa_moves.push_back( StepWiseMove( i_full + 1, Attachment( i_full, BOND_TO_PREVIOUS ), ADD ) );
		}

		// bulge skip...
		if ( ( options_->skip_bulge_frequency() > 0.0 ) &&
				(i == nres || res_list[ i ] + 2 < res_list[ i+1 ]) &&
				i_full < (nres_full - 1)  &&
				!is_cutpoint_in_full_pose[ i_full + 1 ] ) {
			// for now can only handle single residue additions via skip-bulge ("floating base")
			Size const other_pose_idx_intervening = full_model_info.get_idx_for_other_pose_with_residue( i_full + 1 );
			if ( other_pose_idx_intervening == 0 ) {
				if ( is_addable_res( i_full + 2, pose ) ) {
					swa_moves_skip_bulge.push_back( StepWiseMove( i_full + 2, Attachment( i_full, JUMP_TO_PREV_IN_CHAIN ), ADD ) );
				}
			}
		}
	}

	for ( Size i = 1; i <= nres; i++ ) {

		if ( ( i != 1 ) &&
				( !fold_tree.is_cutpoint( i-1 ) || (res_list[ i ]-1 <= res_list[ i-1 ]) ) ) continue;
		// could be a 5' chain terminus

		Size const i_full = res_list[ i ];
		if ( i_full <= 1 || is_cutpoint_in_full_pose[ i_full - 1 ] ) continue;
		// good, there's still a gap!

		// direct addition of single residue (or potentially, a block of residues starting with this residue)
		if ( is_addable_res( i_full - 1, pose ) ) {
			swa_moves.push_back( StepWiseMove( i_full - 1, Attachment( i_full, BOND_TO_NEXT), ADD ) );
		}

		// bulge skip...
		if ( ( options_->skip_bulge_frequency() > 0.0 ) &&
				(i == 1 || res_list[ i ] - 2 > res_list[ i - 1 ]) &&
				i_full > 2 &&
				!is_cutpoint_in_full_pose[ i_full - 2 ] ) {
			// for now can only handle single residue additions via skip-bulge ("floating base")
			Size const other_pose_idx_intervening = full_model_info.get_idx_for_other_pose_with_residue( i_full - 1 );
			if ( other_pose_idx_intervening == 0 ) {
				if ( is_addable_res( i_full - 2, pose ) ) {
					swa_moves_skip_bulge.push_back( StepWiseMove( i_full - 2, Attachment( i_full, JUMP_TO_NEXT_IN_CHAIN ), ADD ) );
				}
			}
		}
	}

	for ( Size n = 1; n <= swa_moves_skip_bulge.size(); n++ ) swa_moves.push_back( swa_moves_skip_bulge[ n ] );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_intramolecular_delete_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	utility::vector1< StepWiseMove > swa_moves_delete;
	get_intramolecular_split_move_elements( pose, swa_moves_delete, DELETE );

	// don't delete the only pose left.
	if ( const_full_model_info( pose ).other_pose_list().size() == 0 && pose.size() == 2 ) return;

	for ( Size n = 1; n <= swa_moves_delete.size(); n++ ) swa_moves.push_back( swa_moves_delete[ n ] );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_intramolecular_split_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves,
	MoveType const move_type ) {

	using namespace protocols::stepwise::modeler::rna;

	utility::vector1< StepWiseMove > swa_moves_split;
	utility::vector1< bool > partition_definition;
	utility::vector1< Size > partition_res1, partition_res2;
	utility::vector1< Size > const domain_map = core::pose::full_model_info::get_input_domain_from_full_model_info_const( pose );
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const chains = figure_out_chain_numbers_from_full_model_info_const( pose );

	// first look at suites
	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( !pose.fold_tree().is_cutpoint( i ) &&
				( domain_map[ i ] != domain_map[ i+1 ] ||
				domain_map[ i ] == 0 || domain_map [ i+1 ] == 0 ) ) {

			partition_definition = get_partition_definition( pose, i /*suite*/ );
			if ( partition_splits_an_input_domain( partition_definition, domain_map ) ) continue;

			partition_res1 = get_partition_res( partition_definition, true );
			partition_res2 = get_partition_res( partition_definition, false );

			// from-scratch handled elsewhere
			if ( move_type == DELETE && both_remnants_would_be_deleted( pose, partition_res1, partition_res2 ) ) continue;
			if ( !options_->allow_submotif_split() && partitions_split_a_submotif( pose, partition_res1, partition_res2 ) ) continue;

			swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full(partition_res2),
				Attachment( full_model_info.sub_to_full(i), BOND_TO_PREVIOUS ), move_type ) );

			if ( !force_unique_moves_ ) { // to match legacy code.
				swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full(partition_res1),
					Attachment( full_model_info.sub_to_full(i+1), BOND_TO_NEXT ), move_type ) );
			}

		}
	}

	// now look at jumps. Note that, currently, additions by jump are only for single residues, so only
	// split cases in which single residues are deleted.
	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {

		partition_definition = get_partition_definition_by_jump( pose, n );
		if ( partition_splits_an_input_domain( partition_definition, domain_map ) ) continue;

		Size const downstream_res = pose.fold_tree().downstream_jump_residue( n );
		Size const upstream_res = pose.fold_tree().upstream_jump_residue( n );

		// by convention, moving res should have higher
		Size const moving_res = std::max( downstream_res, upstream_res );
		Size const reference_res = std::min( downstream_res, upstream_res );
		Size const offset = ( moving_res - reference_res );

		// do not include jumps greater than a single residue (may exist in submotifs)
		if ( std::abs( static_cast< int >( offset ) ) > 2 ) continue;
		if ( chains[ moving_res ] != chains[ reference_res ] ) continue;

		partition_res1 = get_partition_res( partition_definition, ( moving_res < reference_res ) );
		partition_res2 = get_partition_res( partition_definition, ( moving_res > reference_res ) );

		// make sure that at least one partition is a single nucleotide.
		if ( partition_res1.size() > 1 && partition_res2.size() > 1 ) continue;
		if ( !options_->allow_submotif_split() && partitions_split_a_submotif( pose, partition_res1, partition_res2 ) ) continue;

		if ( partition_res1.size() == 1 ) {
			Size const anchor_res = get_anchor_res( partition_res1[1], pose );
			AttachmentType type = ( anchor_res > partition_res1[1] ) ? JUMP_TO_NEXT_IN_CHAIN : JUMP_TO_PREV_IN_CHAIN;
			swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full(partition_res1),
				Attachment( full_model_info.sub_to_full(anchor_res), type ), move_type ) );
			if ( force_unique_moves_ ) continue;
		}
		if ( partition_res2.size() == 1 ) {
			Size const anchor_res = get_anchor_res( partition_res2[1], pose );
			AttachmentType type = ( anchor_res > partition_res2[1] ) ? JUMP_TO_NEXT_IN_CHAIN : JUMP_TO_PREV_IN_CHAIN;
			swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full(partition_res2),
				Attachment( full_model_info.sub_to_full(anchor_res), type ), move_type ) );
		}
	}
	for ( Size n = 1; n <= swa_moves_split.size(); n++ ) swa_moves.push_back( swa_moves_split[ n ] );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::check_from_scratch(  pose::Pose const & pose,
	utility::vector1< bool > const & partition_definition) const {
	if ( options_->from_scratch_frequency() > 0.0 &&
			pose.size() == 2 &&
			!pose.fold_tree().is_cutpoint( 1 ) &&
			partition_definition[ 1 ] != partition_definition[ 2 ] ) {
		return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// docking (cross-chain)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_docking_split_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves,
	MoveType const move_type ){
	using namespace protocols::stepwise::modeler::rna;

	utility::vector1< StepWiseMove > swa_moves_split;
	utility::vector1< bool > partition_definition;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const dock_domain_map = full_model_info.dock_domain_map();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );

	// look at jumps.
	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		partition_definition = get_partition_definition_by_jump( pose, n );
		Size const moving_res = pose.fold_tree().downstream_jump_residue( n );
		Size const reference_res = pose.fold_tree().upstream_jump_residue( n );
		Size const moving_res_full =  res_list[ moving_res ];
		Size const reference_res_full = res_list[ reference_res ];
		if ( !stepwise_addable_pose_residue( moving_res, pose ) ) continue;
		//runtime_assert( dock_domain_map[ moving_res_full ] > 0 );
		//runtime_assert( dock_domain_map[ reference_res_full ] > 0 );
		if ( dock_domain_map[ moving_res_full ] == dock_domain_map[ reference_res_full ] ) continue;
		utility::vector1< Size > const moving_partition_res = get_partition_res( partition_definition,  partition_definition[ moving_res ] );
		utility::vector1< Size > const reference_partition_res = get_partition_res( partition_definition,  partition_definition[ reference_res ] );
		if ( move_type == DELETE && both_remnants_would_be_deleted( pose, moving_partition_res, reference_partition_res ) ) continue;
		if ( !options_->allow_submotif_split() && partitions_split_a_submotif( pose, moving_partition_res, reference_partition_res ) ) continue;

		// order 'moving'/'attachment' based on order in sequence -- help enforce unique StepWiseMove definition.
		if ( moving_res > reference_res ) {
			swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full( moving_partition_res ),
				Attachment( reference_res_full, JUMP_DOCK ), move_type ) );
		} else {
			swa_moves_split.push_back( StepWiseMove( full_model_info.sub_to_full( reference_partition_res ),
				Attachment( moving_res_full, JUMP_DOCK ), move_type ) );
		}
	}

	for ( Size n = 1; n <= swa_moves_split.size(); n++ ) swa_moves.push_back( swa_moves_split[ n ] );
}


//////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::partitions_split_a_submotif(
	pose::Pose const & pose,
	utility::vector1 < Size > const & partition1,
	utility::vector1 < Size > const & partition2
) const {

	utility::vector1< Size > const & pose_res_list( get_res_list_from_full_model_info_const( pose ) );
	utility::vector1< SubMotifInfoOP > submotif_info_list = const_full_model_info( pose ).submotif_info_list();
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list.begin(); itr != submotif_info_list.end(); ++itr ) {
		bool in_partition1( false ), in_partition2( false );
		for ( Size n = 1; n <= (*itr)->res_list().size(); n++ ) {
			Size const submotif_res_in_pose_numbering( pose_res_list.index( (*itr)->res_list( n ) ) );
			if ( partition1.has_value( submotif_res_in_pose_numbering ) ) {
				in_partition1 = true;
				continue;
			}
			if ( partition2.has_value( submotif_res_in_pose_numbering ) ) {
				in_partition2 = true;
				continue;
			}
		}
		if ( in_partition1 && in_partition2 ) {
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::remnant_would_be_deleted(
	pose::Pose const & pose,
	utility::vector1 < Size > const & partition ) const
{
	utility::vector1< Size > const partition_in_full_model_numbering = sub_to_full( partition, pose );
	if ( const_full_model_info( pose ).is_a_submotif( partition_in_full_model_numbering ) &&
			!const_full_model_info( pose ).is_a_submotif_seed( partition_in_full_model_numbering ) ) return true;
	if ( partition.size() == 1 ) return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::both_remnants_would_be_deleted(
	pose::Pose const & pose,
	utility::vector1 < Size > const & partition1,
	utility::vector1 < Size > const & partition2 ) const
{
	return ( remnant_would_be_deleted( pose, partition1 ) && remnant_would_be_deleted( pose, partition2 ) );
}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::check_for_intramolecular_submotif_jump(
	pose::Pose const & pose,
	Size const & moving_res,
	Size const & attached_res
) const {
	utility::vector1< Size > const & check_res = utility::tools::make_vector1( moving_res, attached_res );
	bool const in_a_submotif = const_full_model_info( pose ).in_a_submotif( check_res );
	bool const jump_exists = pose.fold_tree().jump_exists( moving_res, attached_res );
	return ( in_a_submotif && jump_exists );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The point of this function is to allow big move_elements to remain fixed. It may not be entirely necessary.
void
StepWiseMoveSelector::remove_from_consideration_first_multi_residue_move_element( utility::vector1< StepWiseMove > & swa_moves,
	bool remove_even_if_not_singlet ){

	Size number_of_multi_residue_move_elements( 0 ), first_multi_residue_move_element( 0 );
	for ( Size n = 1; n <= swa_moves.size(); n++ ) {
		if ( swa_moves[ n ].move_element().size() > 1 ) {
			number_of_multi_residue_move_elements++;
			if ( first_multi_residue_move_element == 0 ) first_multi_residue_move_element = n;
		}
	}
	if ( number_of_multi_residue_move_elements == 0 ) return;

	runtime_assert( number_of_multi_residue_move_elements > 0 );

	if ( remove_even_if_not_singlet || (number_of_multi_residue_move_elements == 1) ) {
		utility::vector1< StepWiseMove > swa_moves_new;

		for ( Size n = 1; n <= swa_moves.size(); n++ ) {
			if ( n == first_multi_residue_move_element ) continue;
			swa_moves_new.push_back( swa_moves[ n ] );
		}
		swa_moves = swa_moves_new;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// from scratch [create dinucleotides 'anew']
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_from_scratch_add_and_delete_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	get_from_scratch_add_move_elements( pose, swa_moves );
	if ( allow_delete_ ) get_from_scratch_delete_move_elements( pose, swa_moves );
	if ( options_->filter_complex_cycles() ) filter_complex_cycles( swa_moves, pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_from_scratch_add_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ){
	// for RNA, these are dinucleotides that would be 'freely floating' -- not attached
	// to the current pose (or its sisters).
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	Size nres_full = full_model_info.full_sequence().size();

	if ( options_->from_scratch_frequency() == 0.0 ) return;

	Attachments const blank_attachments;

	for ( Size n = 1; n < nres_full; n++ ) {
		if ( cutpoint_open_in_full_model.has_value( n ) ) continue; // must be contiguous dinucleotide.
		if ( already_instantiated_in_pose( pose, n     ) ) continue;
		if ( already_instantiated_in_pose( pose, n + 1 ) ) continue;
		if ( !is_addable_res( n  , pose) ) continue;
		if ( !is_addable_res( n+1, pose ) ) continue;
		swa_moves.push_back( StepWiseMove( utility::tools::make_vector1( n, n+1 ), blank_attachments, FROM_SCRATCH ) );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::figure_out_from_scratch_delete_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) const {

	utility::vector1< Size > const domain_map = core::pose::full_model_info::get_input_domain_from_full_model_info_const( pose );
	utility::vector1< Size > const partition_definition = utility::tools::make_vector1(true,false);
	if ( check_from_scratch( pose, partition_definition ) &&
			!partition_splits_an_input_domain( partition_definition, domain_map ) ) {
		swa_moves.push_back( StepWiseMove( sub_to_full( 2, pose ), Attachment( sub_to_full( 1, pose ), BOND_TO_PREVIOUS ), DELETE ) );
		if ( !force_unique_moves_ ) {
			swa_moves.push_back( StepWiseMove( sub_to_full( 1, pose ), Attachment( sub_to_full( 2, pose ), BOND_TO_NEXT ), DELETE ) );
		}
	}
	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size k = 1; k <= other_pose_list.size(); k++ ) {
		figure_out_from_scratch_delete_move_elements( *( other_pose_list[k] ), swa_moves );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_from_scratch_delete_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ){
	utility::vector1< StepWiseMove > swa_moves_delete;
	figure_out_from_scratch_delete_move_elements( pose, swa_moves_delete );

	// don't delete the only pose left.
	if ( const_full_model_info( pose ).other_pose_list().size() == 0 && pose.size() == 2 ) return;

	for ( Size n = 1; n <= swa_moves_delete.size(); n++ ) swa_moves.push_back( swa_moves_delete[n] );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_docking_add_and_delete_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	get_docking_add_move_elements( pose, swa_moves );
	if ( allow_delete_ ) get_docking_delete_move_elements( pose, swa_moves );
	if ( force_unique_moves_ ) filter_pose_order( swa_moves, pose );
	if ( options_->filter_complex_cycles() ) filter_complex_cycles( swa_moves, pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::figure_out_already_docked(
	pose::Pose const & pose,
	std::map< std::pair< Size, Size >, bool > & already_docked ) const
{
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	std::set< Size > dock_domains = get_dock_domains( res_list, pose );
	for ( auto it1 = dock_domains.begin(),
			end1 = dock_domains.end(); it1 != end1; ++it1 ) {
		for ( auto it2 = it1,
				end2 = dock_domains.end(); it2 != end2; ++it2 ) {
			if ( it1 == it2 ) continue;
			already_docked[ std::make_pair( *it1, *it2 ) ] = true;
			already_docked[ std::make_pair( *it2, *it1 ) ] = true;
		}
	}

	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size k = 1; k <= other_pose_list.size(); k++ ) {
		figure_out_already_docked( *other_pose_list[ k ], already_docked );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::set< Size >
StepWiseMoveSelector::get_dock_domains(
	utility::vector1< Size > const & move_element,
	pose::Pose const & pose ) const
{
	std::set< Size > dock_domains;
	utility::vector1< Size > const & dock_domain_map = const_full_model_info( pose ).dock_domain_map();
	for ( Size i = 1; i <= move_element.size(); i++ ) {
		if ( dock_domain_map[ move_element[ i ] ] > 0 ) dock_domains.insert( dock_domain_map[ move_element[ i ] ] );
	}
	return dock_domains;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_docking_add_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ){

	using namespace core::pose::full_model_info;

	Size const nres( pose.size() );
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	utility::vector1< Size > const & dock_domain_map = const_full_model_info( pose ).dock_domain_map();
	//utility::vector1< Size > const & working_res = const_full_model_info( pose ).working_res();

	if ( res_list.size() == 0 ) return;

	std::map< std::pair< Size, Size >, std::pair< Size, Size > > preferred_jump_pair = get_preferred_jump_pair_for_docking_domains( full_model_info );

	// check for domains that are already docked.
	std::map< std::pair< Size, Size >, bool > already_docked;
	figure_out_already_docked( pose, already_docked );

	for ( Size i = 1; i <= nres; i++ ) {
		Size const i_full = res_list[ i ] ;
		Size const dock_domain_i = dock_domain_map[ i_full ];
		if ( !is_addable_res( i_full, pose ) ) continue;
		for ( Size j_full = 1; j_full <= full_model_info.size(); j_full++ ) {
			if ( res_list.has_value( j_full ) ) continue; // must be docked to different pose.

			Size const dock_domain_j = dock_domain_map[ j_full ];

			if ( already_docked[ std::make_pair( dock_domain_i, dock_domain_j ) ] ) continue;
			if ( dock_domain_j == 0 ) continue; // not a working res.
			if ( dock_domain_i == dock_domain_j ) continue;
			if ( !is_addable_res( j_full, pose ) ) continue;

			if ( preferred_jump_pair.find( std::make_pair( dock_domain_i, dock_domain_j ) ) !=
					preferred_jump_pair.end() ) {
				std::pair< Size, Size > preferred_jump_for_dock_domain = preferred_jump_pair[ std::make_pair( dock_domain_i, dock_domain_j ) ];
				if ( std::make_pair( i_full, j_full ) != preferred_jump_for_dock_domain ) continue;
			}
			swa_moves.push_back( StepWiseMove( j_full, Attachment( i_full, JUMP_DOCK ), ADD ) );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// quick hack -- actually this selects for submotifs that have no intervening bulges *and*
// for which each residue connects to the rest of the pose through some bond -- force base pairs!
bool
check_for_connections_for_all_submotif_residues(
	MoveElement const & move_element,
	utility::vector1< StepWiseMove > const & add_moves )
{
	std::map< Size, bool > addable;
	for ( Size k = 1; k <= add_moves.size(); k++ ) {
		StepWiseMove const & add_move = add_moves[ k ];
		Size const add_move_res = add_move.moving_res();
		addable[ add_move_res ] = true;
	}
	bool all_addable( true );
	for ( Size k = 1; k <= move_element.size(); k++ ) {
		if ( !addable[ move_element[ k ] ] ) all_addable = false;
	}
	return all_addable;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_submotif_add_moves( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ) {

	if ( submotif_library_ == nullptr ) return;

	utility::vector1< StepWiseMove > submotif_from_scratch_moves = submotif_library_->get_submotif_moves( pose );

	// could any of these new submotif_elements by attached to the current pose by a bond?
	// already had developed machinery to figure that out for normal ADD moves:
	utility::vector1< StepWiseMove > add_moves, submotif_add_moves;
	get_intramolecular_add_move_elements( pose, add_moves );
	if ( force_unique_moves_ ) filter_pose_order( add_moves, pose );

	for ( Size n = 1; n <= submotif_from_scratch_moves.size(); n++ ) {
		MoveElement const & move_element = submotif_from_scratch_moves[ n ].move_element();
		std::string const & submotif_tag = submotif_from_scratch_moves[ n ].submotif_tag();

		// this is a hack to force addition of base pair submotifs without intervening bulges.
		if ( options_->force_submotif_without_intervening_bulge() &&
				!check_for_connections_for_all_submotif_residues( move_element, add_moves ) ) continue;

		for ( Size k = 1; k <= add_moves.size(); k++ ) {
			StepWiseMove const & add_move = add_moves[ k ];
			Size const add_move_res = add_move.moving_res();
			if ( move_element.has_value( add_move_res ) ) {
				submotif_add_moves.push_back( StepWiseMove( move_element, add_move.attachments(), ADD_SUBMOTIF, submotif_tag ) );
			}
		} // add_moves
	} // submotifs
	filter_add_submotif_moves_to_not_redock_domain( submotif_add_moves, pose );
	if ( options_->filter_complex_cycles() ) filter_complex_cycles( submotif_add_moves, pose );

	for ( Size n = 1; n <= submotif_add_moves.size(); n++ ) swa_moves.push_back( submotif_add_moves[ n ] );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// don't introduce a second connection between two 'domains' if  there already is one, either via
//  a ADD JUMP_DOCK move or by a previous ADD_SUBMOTIF move.
void
StepWiseMoveSelector::filter_add_submotif_moves_to_not_redock_domain(
	utility::vector1< StepWiseMove > & submotif_add_moves,
	pose::Pose const & pose ) const
{
	utility::vector1< StepWiseMove > swa_moves_new;

	// check for domains that are already docked
	std::map< std::pair< Size, Size >, bool > already_docked;
	figure_out_already_docked( pose, already_docked );

	for ( Size n = 1; n <= submotif_add_moves.size(); n++ ) {
		runtime_assert( submotif_add_moves[ n ].move_type() == ADD_SUBMOTIF );
		// what domains might be docked in submotif?
		std::set< Size > const dock_domains = get_dock_domains( submotif_add_moves[ n ].move_element(), pose );
		// if any are already docked in pose, better not make a second jump connection via this submotif.
		bool ok_to_dock( true );
		for ( auto it1 = dock_domains.begin(),
				end1 = dock_domains.end(); it1 != end1; ++it1 ) {
			for ( auto it2 = it1,
					end2 = dock_domains.end(); it2 != end2; ++it2 ) {
				if ( already_docked[ std::make_pair( *it1, *it2 ) ] ) {
					ok_to_dock = false;
					break;
				}
			}
			if ( !ok_to_dock ) break;
		}
		if ( ok_to_dock ) {
			swa_moves_new.push_back( submotif_add_moves[ n ] );
		}
	}
	submotif_add_moves = swa_moves_new;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
StepWiseMoveSelector::get_unique_chains( pose::Pose const & pose ){
	utility::vector1< Size > const chains = figure_out_chain_numbers_from_full_model_info_const( pose );
	utility::vector1< Size > unique_chains;
	for ( Size n = 1; n <= chains.size(); n++ ) {
		if ( !unique_chains.has_value( chains[n] ) ) unique_chains.push_back( chains[n] );
	}
	return unique_chains;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// this is pretty slow, actually -- can be sped up by pre-computing chains_full, unique chains in each pose, etc.
bool
StepWiseMoveSelector::share_chains( utility::vector1< Size > const & current_unique_chains,
	pose::Pose const & pose,
	Size const j_full ) {

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const chains_full = get_chains_full( pose );

	Size const chain_j = chains_full[ j_full ];
	utility::vector1< Size > other_pose_chains = utility::tools::make_vector1( chain_j );

	Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( j_full );
	if ( other_pose_idx > 0 ) other_pose_chains = get_unique_chains( *full_model_info.other_pose_list()[ other_pose_idx ] );

	for ( Size q = 1; q <= current_unique_chains.size(); q++ ) {
		if ( other_pose_chains.has_value( current_unique_chains[ q ] ) ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_docking_delete_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves ){
	get_docking_split_move_elements( pose, swa_moves, DELETE );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::is_addable_res( Size const n, pose::Pose const & pose ) const {

	// don't allow Mg(2+) or HOH yet -- must be an easier way to figure out ligand or not.
	if ( !stepwise_addable_residue( n, const_full_model_info( pose ).full_model_parameters()->non_standard_residue_map() ) ) return false;
	utility::vector1< Size > const & working_res = const_full_model_info( pose ).working_res();
	utility::vector1< Size > const & bulge_res   = const_full_model_info( pose ).rna_bulge_res();
	return ( working_res.has_value( n ) && !bulge_res.has_value( n ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Note that loops are returned in parent full_model numbering, i.e. numbering assuming maximal
// loop lengths.
// Following assumes that loops cannot be split into multiple loops, e.g., through submotif moves...
// That case could be handled too.
utility::vector1< std::pair< Size, Size > >
get_design_loops( FullModelParametersCOP parent_full_model_parameters )
{
	utility::vector1< std::pair< Size, Size > > design_loops;
	utility::vector1< Size > const & cutpoint_open = parent_full_model_parameters->get_res_list( CUTPOINT_OPEN );
	std::string const & full_sequence( parent_full_model_parameters->full_sequence() );
	bool in_loop( false );
	Size start_loop( 1 );
	Size const nres( full_sequence.size() );
	for ( Size n = 1; n <= nres + 1; n++ ) {
		bool const is_designable_residue = ( n <= nres && full_sequence[ n - 1 ] == 'n' );
		if ( in_loop && ( !is_designable_residue || cutpoint_open.has_value( n - 1 ) ) ) {
			in_loop = false;
			Size const end_loop = n - 1;
			runtime_assert( end_loop >= start_loop );
			design_loops.push_back( std::make_pair( start_loop, end_loop ) );
		}
		if ( !in_loop && is_designable_residue ) {
			in_loop = true;
			start_loop = n;
		}
	}
	return design_loops;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_vary_loop_length_moves( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & vary_loop_length_moves ) const
{
	FullModelParametersCOP full_model_parameters( const_full_model_info( pose ).full_model_parameters() );
	FullModelParametersCOP parent_full_model_parameters( full_model_parameters->parent_full_model_parameters() );
	runtime_assert( parent_full_model_parameters != nullptr );
	utility::vector1< Size > const & slice_res_list = const_full_model_info( pose ).full_model_parameters()->slice_res_list();
	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();
	Size const nres = full_model_parameters->size();

	// loops in parent full model numbering:
	utility::vector1< std::pair< Size, Size > > design_loops = get_design_loops( parent_full_model_parameters );
	for ( Size n = 1; n <= design_loops.size(); n++ ) {
		std::pair< Size, Size > const & design_loop = design_loops[ n ];
		Size const max_gap_size = design_loop.second - design_loop.first + 1;
		runtime_assert( max_gap_size > 0 );

		// following are in full_model numbering.
		// takeoff_res is position immediately previous to loop residues, landing_res is immediately after.
		Size takeoff_res( 0 ), landing_res( 0 );
		if ( slice_res_list.has_value( design_loop.first - 1 ) ) {
			takeoff_res = slice_res_list.index( design_loop.first - 1 );
		} else {
			bool found_takeoff_res( false );
			for ( Size q = design_loop.first; q <= design_loop.second + 1; q++ ) {
				if ( slice_res_list.has_value( q ) ) {
					takeoff_res = slice_res_list.index( q ) - 1;
					found_takeoff_res = true;
					break;
				}
			}
			runtime_assert( found_takeoff_res );
		}

		if ( slice_res_list.has_value( design_loop.second + 1 ) ) {
			landing_res = slice_res_list.index( design_loop.second + 1 );
		} else {
			bool found_landing_res( false );
			for ( Size q = design_loop.second; q >= design_loop.first - 1; q-- ) {
				if ( slice_res_list.has_value( q ) ) {
					landing_res = slice_res_list.index( q ) + 1;
					found_landing_res = true;
					break;
				}
			}
			runtime_assert( found_landing_res );
		}

		runtime_assert( takeoff_res < landing_res );
		Size const gap_size = landing_res - takeoff_res - 1;
		runtime_assert( gap_size <= max_gap_size );

		bool const tethered_at_takeoff = ( takeoff_res >= 1    && !cutpoint_open_in_full_model.has_value( takeoff_res ) );
		bool const tethered_at_landing = ( landing_res <= nres && !cutpoint_open_in_full_model.has_value( landing_res - 1 ) );
		runtime_assert( tethered_at_takeoff || tethered_at_landing );

		// can we remove a nucleotide from the loop? Must be at a terminal
		if ( gap_size > 0 ) {
			// can either remove from beginning or end of loop -- move is equivalent. by convention, remove at end, unless
			// beginning is 'free', in which case remove from there.
			if ( tethered_at_takeoff ) {
				vary_loop_length_moves.push_back( StepWiseMove( landing_res-1, Attachment( landing_res-2, BOND_TO_PREVIOUS ), DELETE_LOOP_RES ) );
			} else {
				vary_loop_length_moves.push_back( StepWiseMove( takeoff_res+1, Attachment( takeoff_res+2, BOND_TO_NEXT ), DELETE_LOOP_RES ) );
			}
		}

		// can we insert a nucleotide into the loop? Again, must be inserted at a terminal
		if ( gap_size < max_gap_size ) {
			// can either insert into beginning or end of loop -- move is equivalent. by convention, insert at end, unless
			// beginning is 'free', in which case insert there.
			if ( tethered_at_takeoff ) {
				vary_loop_length_moves.push_back( StepWiseMove( landing_res, Attachment( landing_res-1, BOND_TO_PREVIOUS ), ADD_LOOP_RES ) );
			} else {
				vary_loop_length_moves.push_back( StepWiseMove( takeoff_res, Attachment( takeoff_res+1, BOND_TO_NEXT ), ADD_LOOP_RES ) );
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// don't break up any fixed domains by a split.
bool
StepWiseMoveSelector::partition_splits_an_input_domain( utility::vector1< Size > const & partition_definition,
	utility::vector1< Size > const & domain_map ) const {
	// parse out which residues go into each domain.
	utility::vector1< Size > blank_vector;
	utility::vector1< utility::vector1< Size > > all_domain_res;
	for ( Size i = 1; i <= domain_map.size(); i++ ) {
		Size const & domain = domain_map[ i ];
		for ( Size n = all_domain_res.size(); n < domain; n++ ) all_domain_res.push_back( blank_vector );
		if ( domain > 0 ) all_domain_res[ domain ].push_back( i );
	}
	for ( Size n = 1; n <= all_domain_res.size(); n++ ) {
		// for this domain, check if all residues are in a single partition.
		bool in_partition_0( false ), in_partition_1( false );
		for ( Size k = 1; k <= all_domain_res[ n ].size(); k++ ) {
			if ( partition_definition[ all_domain_res[ n ][ k ] ] ) {
				in_partition_0 = true;
			} else {
				in_partition_1 = true;
			}
			if ( in_partition_0 && in_partition_1 ) return true;
		}
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_terminal_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves,
	MoveType const & move_type ) {

	using namespace core::pose::full_model_info;

	if ( get_res_list_from_full_model_info_const( pose ).size() == 0 ) return;
	utility::vector1< MoveElement > const & move_elements = get_move_elements_from_full_model_info_const( pose );

	Size const num_move_elements( move_elements.size() );
	for ( Size n = 1; n <= num_move_elements; n++ ) {

		MoveElement const & move_element = move_elements[ n ];
		utility::vector1< Attachment > attachments = get_attachments( pose, move_element );
		// at least for now, each move_element should be connected to another one.
		if ( attachments.empty() ) {
			runtime_assert( num_move_elements == 1 );
			continue;
		}
		if ( attachments.size() > 1 ) continue; // not a terminal
		// for now, floating base can only handle single residues.
		if ( attachments.size() == 1 &&
				( attachments[1].attachment_type() == JUMP_TO_PREV_IN_CHAIN || attachments[1].attachment_type() == JUMP_TO_NEXT_IN_CHAIN ) &&
				move_element.size() > 1 ) continue;

		swa_moves.emplace_back( move_element, attachments, move_type );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::get_internal_move_elements( pose::Pose const & pose,
	utility::vector1< StepWiseMove > & swa_moves,
	MoveType const & move_type ) {

	using namespace core::pose::full_model_info;

	utility::vector1< MoveElement > const & move_elements = get_move_elements_from_full_model_info_const( pose );

	Size const num_move_elements( move_elements.size() );
	for ( Size n = 1; n <= num_move_elements; n++ ) {

		MoveElement const & move_element = move_elements[ n ];
		utility::vector1< Attachment > attachments = get_attachments( pose, move_element );

		// at least for now, each move_element should be connected to another one.
		if ( attachments.size() != 2 ) continue;

		// note that PREVIOUS should be before NEXT.
		if ( attachments[1].attachment_type() == BOND_TO_PREVIOUS &&
				attachments[2].attachment_type() == BOND_TO_NEXT ) {
			swa_moves.push_back( StepWiseMove( move_element, attachments, move_type ) );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
StepWiseMoveSelector::get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting ) const {
	utility::vector1< Size > partition_res;
	for ( Size n = 1; n <= partition_definition.size(); n++ ) {
		if ( partition_definition[ n ]  == setting ) partition_res.push_back( n );
	}
	return partition_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Attachments
StepWiseMoveSelector::get_attachments( pose::Pose const & pose, Size const & moving_res ){
	return get_attachments( pose, utility::tools::make_vector1( moving_res ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A moving element (or even a single residue ) can actually have several connections.
// Some are covalent (polymeric) connections (BOND_TO_PREVIOUS, BOND_TO_NEXT)
// Some could involve jumps to residues in same chain (e.g., in 'skip-bulge' moves).
// Some could involve jumps to other chains or ligands.
Attachments
StepWiseMoveSelector::get_attachments( pose::Pose const & pose, MoveElement const & move_element ){

	using namespace core::pose::full_model_info;

	Size const & nres( pose.size() );
	kinematics::FoldTree const & fold_tree( pose.fold_tree() );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );

	Attachments attachments_to_previous, attachments_to_next, attachments;

	// need to find number of attachment points of this move_element to rest of pose.
	for ( Size k = 1; k <= move_element.size(); k++ ) {

		Size const i_full = move_element[ k ];
		Size const i = res_list.index( i_full );

		if ( i > 1 && ! move_element.has_value( res_list[ i - 1 ] ) ) {
			// potential connection to previous residue in rest of pose
			if ( !fold_tree.is_cutpoint( i - 1 ) )  {
				attachments_to_previous.push_back( Attachment( i_full - 1, BOND_TO_PREVIOUS ) );
			}
			// potential connection to previous residue by jump (a result of 'skip bulge' move)
			Size const jump_attachment_res = check_jump_to_previous_residue_in_chain( pose, i, move_element, const_full_model_info( pose ) );
			if ( jump_attachment_res > 0 ) attachments_to_previous.push_back( Attachment( res_list[ jump_attachment_res ], JUMP_TO_PREV_IN_CHAIN ) );
		}

		// potential connection to next residue in rest of pose
		if ( i < nres && ! move_element.has_value( res_list[ i + 1 ] ) ) {
			if ( !fold_tree.is_cutpoint( i ) ) {
				attachments_to_next.push_back( Attachment( i_full + 1, BOND_TO_NEXT ) );
			}
			// potential connection to next by jump (a result of 'skip bulge' move)
			Size const jump_attachment_res = check_jump_to_next_residue_in_chain( pose, i, move_element, const_full_model_info( pose ) );
			if ( jump_attachment_res > 0 ) attachments_to_next.push_back( Attachment( res_list[ jump_attachment_res ], JUMP_TO_NEXT_IN_CHAIN ) );
		}
	}

	for ( Size k = 1; k <= attachments_to_previous.size(); k++ ) attachments.push_back( attachments_to_previous[ k ] );
	for ( Size k = 1; k <= attachments_to_next.size()    ; k++ ) attachments.push_back( attachments_to_next[ k ] );

	// have to relax this when we enable inter-chain jumps, or enable multiple intra-chain jumps.
	//  runtime_assert( attachments.size() <= 4 );
	return attachments;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMoveSelector::already_instantiated_in_pose( pose::Pose const & pose, Size const & resnum_in_full_model_numbering ) const {

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	if ( res_list.size() == 0 ) return false;

	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( res_list[ n ] == resnum_in_full_model_numbering ) return true;
	}

	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size k = 1; k <= other_pose_list.size(); k++ ) {
		if ( already_instantiated_in_pose( *(other_pose_list[ k ]), resnum_in_full_model_numbering ) ) return true;
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
StepWiseMoveSelector::figure_out_first_res_in_pose( utility::vector1< Size > const & pose_domain_map ) const {
	utility::vector1< Size > first_res_in_pose( max( pose_domain_map ), 0 );
	for ( Size n = 1; n <= pose_domain_map.size(); n++ ) {
		if ( pose_domain_map[ n ] > 0  && first_res_in_pose[ pose_domain_map[ n ] ] == 0 ) {
			first_res_in_pose[ pose_domain_map[ n ] ] = n;
		}
	}
	return first_res_in_pose;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::filter_pose_order( utility::vector1< StepWiseMove > & swa_moves,
	pose::Pose const & pose ) const {

	runtime_assert( force_unique_moves_ );

	utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map_const( pose );
	utility::vector1< Size > const first_res_in_pose = figure_out_first_res_in_pose( pose_domain_map );
	utility::vector1< Size > const & sample_res = const_full_model_info( pose ).sample_res();

	utility::vector1< StepWiseMove > swa_moves_filtered;
	for ( Size i = 1; i <= swa_moves.size(); i++ ) {

		StepWiseMove const & swa_move = swa_moves[ i ];

		if ( swa_move.move_type() == ADD ) {
			Size const moving_res   = swa_move.moving_res();
			Size const attached_res = swa_move.attached_res();

			Size const pose_idx1 = pose_domain_map[ moving_res ];
			Size const pose_idx2 = pose_domain_map[ attached_res ];

			if ( pose_idx1 > 0 && pose_idx2 > 0 &&
					pose_idx1 != pose_idx2 &&
					sample_res.has_value( attached_res ) && // check that reverse move is really viable
					first_res_in_pose[ pose_idx1 ] > first_res_in_pose[ pose_idx2 ] ) continue;
		}
		swa_moves_filtered.push_back( swa_move );
	}

	swa_moves = swa_moves_filtered;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Need to do some serious cleanup after getting this to compile...
StepWiseMove
StepWiseMoveSelector::reverse_move( StepWiseMove const & swa_move, pose::Pose const & pose_before, pose::Pose const & pose_after ) const {
	MoveType const & move_type = swa_move.move_type();

	if ( move_type == ADD ) {
		return reverse_add_move( swa_move, pose_after );
	} else if ( move_type == DELETE ) {
		return reverse_delete_move( swa_move, pose_before, pose_after);
	} else if ( move_type == FROM_SCRATCH ) {
		return StepWiseMove( swa_move.move_element()[2], Attachment( swa_move.move_element()[1], BOND_TO_PREVIOUS ), DELETE );
	} else if ( move_type == RESAMPLE ) {
		return reverse_resample_move( swa_move, pose_after );
	} else if ( move_type == RESAMPLE_INTERNAL_LOCAL ) {
		return swa_move;
	} else if ( move_type == ADD_SUBMOTIF ) {
		return reverse_add_submotif_move( swa_move, pose_after );
	} else if ( move_type == ADD_LOOP_RES ) {}
	else if ( move_type == DELETE_LOOP_RES ) {}
	else {
		runtime_assert( move_type == NO_MOVE );
	}
	return StepWiseMove(); // empty
}

//////////////////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseMoveSelector::get_actual_moving_res( StepWiseMove const & swa_move, pose::Pose const & pose ) const {

	FullModelInfo const & full_model_info = const_full_model_info( pose );

	// need to watch out for any kind of rerooting or jump sliding.
	// first break into partitions.
	utility::vector1< Size >  moving_partition_res = full_model_info.full_to_sub( swa_move.move_element() );
	utility::vector1< Size >  root_partition_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !moving_partition_res.has_value( n ) ) root_partition_res.push_back( n );
	}

	// rerooting may have occurred between poses.
	Size actual_moving_res = find_downstream_connection_res( pose, moving_partition_res );
	if ( actual_moving_res == 0 ) actual_moving_res = find_downstream_connection_res( pose, root_partition_res );
	return actual_moving_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseMove
StepWiseMoveSelector::reverse_delete_move( StepWiseMove const & swa_move, pose::Pose const & pose_before, pose::Pose const & pose_after ) const {

	Size pose_res = swa_move.move_element()[ 1 ];
	Size const other_pose_idx = const_full_model_info( pose_before ).get_idx_for_other_pose_with_residue( pose_res );
	Pose const & old_pose = ( other_pose_idx == 0 ) ? pose_before :
		*const_full_model_info( pose_before ).other_pose_list()[ other_pose_idx ];

	// pair of residues that were actually connected in old_pose before the split.
	Size const actual_moving_res = get_actual_moving_res( swa_move, old_pose ); // in working numbering
	Size const res1 = sub_to_full( actual_moving_res, old_pose );
	Size const res2 = sub_to_full( old_pose.fold_tree().get_parent_residue( actual_moving_res ), old_pose );

	utility::vector1< Size > pose_domain_map = figure_out_pose_domain_map_const( pose_after );
	Size const & pose_domain1 = pose_domain_map[ res1 ];
	Size const & pose_domain2 = pose_domain_map[ res2 ];
	if ( pose_domain1 > 0 && pose_domain2 > 0 ) {
		utility::vector1< Size > const first_res_in_pose = figure_out_first_res_in_pose( pose_domain_map );
		utility::vector1< Size > const & sample_res = const_full_model_info( pose_before ).sample_res();
		if ( ( first_res_in_pose[ pose_domain1 ] < first_res_in_pose[ pose_domain2 ] && sample_res.has_value( res1 ) ) ||
				!sample_res.has_value( res2 ) /*special case -- adding a chunk to a fixed input pose*/ ) {
			return StepWiseMove( res1, figure_out_attachment( res1, res2, old_pose ), ADD );
		} else {
			return StepWiseMove( res2, figure_out_attachment( res2, res1, old_pose ), ADD );
		}
	} else if ( pose_domain1 == 0 && pose_domain2 > 0 ) {
		return StepWiseMove( res1, figure_out_attachment( res1, res2, old_pose ), ADD );
	} else if ( pose_domain2 == 0 && pose_domain1 > 0 ) {
		return StepWiseMove( res2, figure_out_attachment( res2, res1, old_pose ), ADD );
	} else {
		runtime_assert( pose_domain1 == 0 && pose_domain2 == 0 );
		if ( old_pose.size() != 2 || !old_pose.fold_tree().is_cutpoint( 1 ) ) {
			TR << "Was expecting old_pose to be a result of from_scratch, with two residues and no cutpoints" << std::endl;
			TR << "old_pose res_list: " << const_full_model_info( old_pose ).res_list() << "  and fold_tree " << old_pose.fold_tree() << std::endl;
			TR << "old_pose submotif_info_list: ";
			const_full_model_info( old_pose ).show_submotif_info_list();
		}
		runtime_assert( old_pose.size() == 2 && !old_pose.fold_tree().is_cutpoint( 1 ) );
		runtime_assert( std::abs( int( res1 ) - int( res2 ) ) == 1 );
		return StepWiseMove( utility::tools::make_vector1( std::min(res1,res2), std::max(res1,res2) ), Attachments(), FROM_SCRATCH );
	}

	//fpd  dummy statement for llvm-gcc-4.2 warnings-as-errors
	return StepWiseMove();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseMove
StepWiseMoveSelector::reverse_add_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const {

	FullModelInfo const & full_model_info = const_full_model_info( pose_after );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	runtime_assert( res_list.has_value( swa_move.moving_res() ) && res_list.has_value( swa_move.attached_res() ) );

	Size working_moving_res   = res_list.index( swa_move.moving_res() );
	Size working_attached_res = res_list.index( swa_move.attached_res() );
	Size actual_moving_res, actual_attached_res;
	if ( Size( pose_after.fold_tree().get_parent_residue( working_moving_res ) ) ==  working_attached_res ) {
		actual_moving_res = working_moving_res;
		actual_attached_res = working_attached_res;
	} else {
		// may have re-rooted
		if ( Size( pose_after.fold_tree().get_parent_residue( working_attached_res ) ) != working_moving_res )  {
			TR << "Res list in pose (full-model): " << res_list << std::endl;
			TR << "Checking that moving_res " << swa_move.moving_res() << " is parent of attached_res " << swa_move.attached_res() << std::endl;
			TR << "Checking that working_moving_res " << working_moving_res << " is parent of working_attached_res " << working_attached_res << std::endl;
			TR << "Fold tree (working-numbering): " << pose_after.fold_tree() << std::endl;
		}
		runtime_assert( Size( pose_after.fold_tree().get_parent_residue( working_attached_res ) ) == working_moving_res );
		actual_moving_res   = working_attached_res;
		actual_attached_res = working_moving_res;
	}

	return ordered_move_from_partition( actual_moving_res, actual_attached_res, pose_after, DELETE );
}


//////////////////////////////////////////////////////////////////////////////////////////////
// NOTE -- the ADD_SUBMOTIF move is now reversible.  It is really two moves.
//  It (1) creates a submotif anew, and (2) adds the submotif to the existing pose.
// However, we now recorded SubMotifInfo in full_model_info, so that if a submotif is
//  split off, it is also deleted -- providing an exact reverse!
//////////////////////////////////////////////////////////////////////////////////////////////
StepWiseMove
StepWiseMoveSelector::reverse_add_submotif_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const {
	Size const actual_moving_res   = get_actual_moving_res( swa_move, pose_after );
	Size const actual_attached_res = pose_after.fold_tree().get_parent_residue( actual_moving_res );
	return ordered_move_from_partition( actual_moving_res, actual_attached_res, pose_after, DELETE );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseMove
StepWiseMoveSelector::reverse_resample_move( StepWiseMove const & swa_move, pose::Pose const & pose_after ) const {

	Size actual_moving_res = get_actual_moving_res( swa_move, pose_after );
	Size const actual_attached_res = pose_after.fold_tree().get_parent_residue( actual_moving_res );
	return ordered_move_from_partition( actual_moving_res, actual_attached_res, pose_after, RESAMPLE );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseMove
StepWiseMoveSelector::ordered_move_from_partition( Size const actual_moving_res, Size const actual_attached_res,
	pose::Pose const & pose, MoveType const & move_type ) const {

	utility::vector1< Size > root_partition_res, moving_partition_res;
	figure_out_root_and_moving_partition_res( pose, actual_moving_res, root_partition_res, moving_partition_res );

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	MoveElement move_element;
	Size moving_res, attached_res;
	// convention used in force_unique_moves_through_pose_order...
	if ( actual_moving_res > actual_attached_res ) {
		move_element = full_model_info.sub_to_full( moving_partition_res );
		moving_res   = full_model_info.sub_to_full( actual_moving_res );
		attached_res = full_model_info.sub_to_full( actual_attached_res );
	} else {
		move_element = full_model_info.sub_to_full( root_partition_res );
		moving_res   = full_model_info.sub_to_full( actual_attached_res );
		attached_res = full_model_info.sub_to_full( actual_moving_res );
	}

	Attachment attachment = figure_out_attachment( moving_res, attached_res, pose );
	return StepWiseMove( move_element, attachment, move_type );
}


//////////////////////////////////////////////////////////////////////////////////////////////////
Attachment
StepWiseMoveSelector::figure_out_attachment( Size const moving_res, Size const attached_res,
	pose::Pose const & pose ) const {

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	utility::vector1< Size > const chains_full = get_chains_full( pose );
	utility::vector1< Size > const dock_domain_map = full_model_info.dock_domain_map();

	Attachment attachment;
	if ( moving_res == attached_res + 1 && !cutpoint_open_in_full_model.has_value( attached_res ) ) {
		attachment = Attachment( attached_res, BOND_TO_PREVIOUS );
	} else if ( moving_res == attached_res - 1 && !cutpoint_open_in_full_model.has_value( moving_res ) ) {
		attachment = Attachment( attached_res, BOND_TO_NEXT );
	} else if ( ( moving_res == attached_res + 2 &&
			chains_full[ moving_res ] == chains_full[ attached_res ] ) ||
			( moving_res > attached_res + 2 &&
			check_for_intramolecular_submotif_jump( pose, moving_res, attached_res ) ) ) {
		attachment = Attachment( attached_res, JUMP_TO_PREV_IN_CHAIN );
	} else if ( ( moving_res == attached_res - 2 &&
			chains_full[ moving_res ] == chains_full[ attached_res ] ) ||
			( moving_res < attached_res - 2 &&
			check_for_intramolecular_submotif_jump( pose, moving_res, attached_res ) ) ) {
		attachment = Attachment( attached_res, JUMP_TO_NEXT_IN_CHAIN );
	} else {
		//runtime_assert( dock_domain_map[ moving_res ] != dock_domain_map[ attached_res ] );
		attachment = Attachment( attached_res, JUMP_DOCK );
	}
	return attachment;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Currently loop_close score term can only handle models whose uninstantiated pose-to-pose loop
//   interconnections form  'simple' cycles -- no cycles can share a loop. See core/scoring/LoopGraph.cc
//   for an explanation.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMoveSelector::filter_complex_cycles( utility::vector1< StepWiseMove > & swa_moves, pose::Pose const & pose ) const {
	utility::vector1< StepWiseMove > swa_moves_new;
	for ( Size n = 1; n <= swa_moves.size(); n++ ) {
		if ( just_simple_cycles( swa_moves[ n ], pose ) ) swa_moves_new.push_back( swa_moves[ n ] );
	}
	swa_moves = swa_moves_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// DEPRECATE verbose if not in use after March 2015.
bool
StepWiseMoveSelector::just_simple_cycles( StepWiseMove const & swa_move, pose::Pose const & pose, bool const verbose /* = false */ ) const {

	utility::vector1< Size > pose_domain_map = figure_out_pose_domain_map_const( pose );
	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();

	MoveType const & move_type = swa_move.move_type();
	MoveElement const & move_element = swa_move.move_element();

	if ( move_type == ADD ) {
		Size const old_domain = pose_domain_map[ swa_move.moving_res() ];
		Size const new_domain = pose_domain_map[ swa_move.attached_res() ];
		if ( old_domain == 0 ) { // building new residue
			pose_domain_map[ swa_move.moving_res() ] = new_domain;
		} else { // this is a *merge* of two existing poses.
			for ( Size n = 1; n <= pose_domain_map.size(); n++ ) {
				if ( pose_domain_map[ n ] == old_domain ) pose_domain_map[ n ] = new_domain;
			}
		}
	} else if ( move_type == DELETE ) {
		// need to figure out if either partition will disappear.
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > other_element;
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			if ( !move_element.has_value( res_list[ n ] ) ) other_element.push_back( res_list[ n ] );
		}
		if ( const_full_model_info( pose ).is_a_submotif( move_element ) &&
				!const_full_model_info( pose ).is_a_submotif_seed( move_element ) ) {
			for ( Size k = 1; k <= move_element.size(); k++ ) pose_domain_map[ move_element[ k ] ] = 0;
		} else if (  const_full_model_info( pose ).is_a_submotif( other_element ) &&
				!const_full_model_info( pose ).is_a_submotif_seed( other_element ) ) {
			for ( Size k = 1; k <= other_element.size(); k++ ) pose_domain_map[ other_element[ k ] ] = 0;
		} else {
			Size const new_domain = max( pose_domain_map ) + 1;
			for ( Size k = 1; k <= move_element.size(); k++ ) pose_domain_map[ move_element[ k ] ] = new_domain;
		}
	} else if ( move_type == FROM_SCRATCH ) {
		Size const new_domain = max( pose_domain_map ) + 1;
		for ( Size k = 1; k <= move_element.size(); k++ ) pose_domain_map[ move_element[ k ] ] = new_domain;
	} else if ( move_type == RESAMPLE || move_type == RESAMPLE_INTERNAL_LOCAL ) {
		// no op
	} else if ( move_type == ADD_SUBMOTIF ) {
		Size const new_domain = pose_domain_map[ swa_move.attached_res() ];
		for ( Size k = 1; k <= move_element.size(); k++ ) pose_domain_map[ move_element[ k ] ] = new_domain;
	} else {
		runtime_assert( move_type == NO_MOVE );
	}

	// Using this LoopGraph object to check cycles
	// Its the same object actually used in the scorefunction to compute loop_close,
	//  which doesn't work when loop graph has cycles that share loops -- too complex
	//  a network to analytically calculate an entropy loss from the connections.
	// Could have StepWiseMoveSelector hold this object as private data for speed, but
	//  I like that its held in here like a utility.
	core::scoring::loop_graph::LoopGraph loop_graph;
	loop_graph.set_error_out_on_complex_cycles( false );
	loop_graph.update_loops_and_cycles( pose_domain_map, cutpoint_open_in_full_model );
	if ( verbose ) { // deprecate on April 2015 or later if not in use.
		TR << "CHECKING SIMPLE CYCLES " << swa_move << std::endl;
		TR << pose_domain_map << std::endl;
		loop_graph.check_loop_cycles_are_disjoint( true );
	}
	return loop_graph.has_just_simple_cycles();
}

} //mover
} //monte_carlo
} //stepwise
} //protocols
