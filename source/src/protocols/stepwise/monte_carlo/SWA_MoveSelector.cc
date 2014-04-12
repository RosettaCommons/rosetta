// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/SWA_MoveSelector.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh> // probably do not need RNA in here.
#include <protocols/moves/MonteCarlo.hh>
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/RNA_Util.hh> // Probably could get rid of RNA-specialized stuff here.
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

#include <map>
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG(239111);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.stepwise.monte_carlo.SWA_MoveSelector" );

using namespace core;
using namespace core::pose::full_model_info;

namespace protocols {
namespace stepwise {
namespace monte_carlo {

	//Constructor
	SWA_MoveSelector::SWA_MoveSelector():
		allow_delete_( true ),
		allow_skip_bulge_( false ),
		from_scratch_frequency_( 0.0 ),
		intermolecular_frequency_( 0.0 ),
		only_dock_preexisting_chunks_( false ),
		allow_shared_chains_in_dock_poses_( false ),
		resampling_( false )
	{}

	//Destructor
	SWA_MoveSelector::~SWA_MoveSelector()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// The reason to consider add and delete together is to help ensure detailed balance.
	//
	// The current scheme (intermolecular_frequency_ ) allows user input of how often docking vs. folding is sampled.
	//  however it does not quite follow detailed balance. Instead, would need to compile all add & delete
	//  moves, and then compute rates for each move, weighted relative to, e.g., the addition of a single residue.
	//  Not all deletes are the same! [their rates would depend on our rate scheme for additions].
	//  I have worked this out, but will need to implement carefully. -- rhiju, march 2014.
	//
	void
	SWA_MoveSelector::get_random_add_or_delete_element( pose::Pose const & pose,
																											SWA_Move & swa_move,
																											utility::vector1< Size > const & sample_res /*leave empty if no filter*/) {
		using namespace core::pose::full_model_info;

		utility::vector1< SWA_Move > swa_moves;
		if ( RG.uniform() < intermolecular_frequency_ ) {
			get_intermolecular_add_and_delete_elements( pose, swa_moves, sample_res ); // docking
		} else if ( RG.uniform() < from_scratch_frequency_ ){
			get_from_scratch_add_and_delete_elements( pose, swa_moves, sample_res ); // docking
		}
		if ( swa_moves.size() == 0  ) get_intramolecular_add_and_delete_elements( pose, swa_moves, sample_res ); // folding

		// backup -- in case we're starting really from nil.
		if ( swa_moves.size() == 0  ) get_from_scratch_add_and_delete_elements( pose, swa_moves, sample_res ); // docking
		if ( swa_moves.size() == 0  ) get_intermolecular_add_and_delete_elements( pose, swa_moves, sample_res ); // folding

		for ( Size n = 1; n <= swa_moves.size(); n++ )	TR.Debug << TR.Green << swa_moves[ n ] << TR.Reset << std::endl;

		if ( swa_moves.size() == 0 ){
			swa_move = SWA_Move();
			return;
		}

		swa_move = RG.random_element( swa_moves );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_random_add_or_delete_element( pose::Pose const & pose,
																											SWA_Move & swa_move ){

		utility::vector1< Size > sample_res; /*leave empty if no filter*/
		get_random_add_or_delete_element( pose, swa_move,
																			sample_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_resample_move_elements( pose::Pose const & pose,
																								utility::vector1< SWA_Move > & swa_moves ) {
		resampling_ = true;
		if ( allow_internal_hinge_ ){
			get_resample_internal_move_elements( pose, swa_moves );
		} else {
			get_resample_terminal_move_elements( pose, swa_moves );
		}
		for ( Size n = 1; n <= swa_moves.size(); n++ ){
			TR.Debug  << swa_moves[ n ] << std::endl;
		}
		if ( allow_internal_local_ ){
			get_resample_internal_local_move_elements( pose, swa_moves );
		}
		resampling_ = false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_resample_terminal_move_elements( pose::Pose const & pose,
																												 utility::vector1< SWA_Move > & swa_moves ) {

		utility::vector1< SWA_Move > swa_moves_terminal;
		get_terminal_move_elements( pose, swa_moves_terminal, RESAMPLE );
		remove_from_consideration_first_multi_residue_move_element( swa_moves_terminal, true /*remove_even_if_not_singlet*/ );
		for ( Size n = 1; n <= swa_moves_terminal.size(); n++ ) swa_moves.push_back( swa_moves_terminal[ n ] );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_resample_internal_move_elements( pose::Pose const & pose,
																												 utility::vector1< SWA_Move > & swa_moves ) {
		utility::vector1< SWA_Move > swa_moves_internal;
		get_intramolecular_split_move_elements( pose, swa_moves_internal, RESAMPLE );
		if ( intermolecular_frequency_ > 0.0 ) get_intermolecular_split_move_elements( pose, swa_moves_internal, RESAMPLE );
		for ( Size n = 1; n <= swa_moves_internal.size(); n++ ) swa_moves.push_back( swa_moves_internal[ n ] );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_resample_internal_local_move_elements( pose::Pose const & pose,
																															 utility::vector1< SWA_Move > & swa_moves ) {

		utility::vector1< SWA_Move > swa_moves_internal;
		get_internal_move_elements( pose, swa_moves_internal, RESAMPLE_INTERNAL_LOCAL );
		// don't delete a multi_residue_move_element if its the only one!
		remove_from_consideration_first_multi_residue_move_element( swa_moves_internal, false /*remove_even_if_not_singlet*/ );
		for ( Size n = 1; n <= swa_moves_internal.size(); n++ ) swa_moves.push_back( swa_moves_internal[ n ] );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// intramolecular [in-chain]
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	SWA_MoveSelector::get_intramolecular_add_and_delete_elements( pose::Pose const & pose,
																							utility::vector1< SWA_Move > & swa_moves,
																							utility::vector1< Size > const & sample_res /*leave empty if no filter*/) {

		if ( allow_delete_ )     get_intramolecular_delete_move_elements( pose, swa_moves );
		get_intramolecular_add_move_elements( pose, swa_moves );
		if ( sample_res.size() > 0 ) filter_by_sample_res( swa_moves, sample_res );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intramolecular_add_move_elements( pose::Pose const & pose,
																													utility::vector1< SWA_Move > & swa_moves ){

		using namespace core::pose::full_model_info;

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		Size nres_full = full_model_info.full_sequence().size();

		if ( res_list.size() == 0 ) return;

		utility::vector1< bool > is_cutpoint_in_full_pose;
		for ( Size i = 1; i <= nres_full; i++ ) is_cutpoint_in_full_pose.push_back( false );
		for ( Size n = 1; n <= cutpoint_open_in_full_model.size(); n++ ) is_cutpoint_in_full_pose[ cutpoint_open_in_full_model[n] ] = true;

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == nres ) ||
					 ( fold_tree.is_cutpoint( i ) && (res_list[ i ]+1 < res_list[ i+1 ]) ) ) { // could be a 3' chain terminus

				Size const i_full = res_list[ i ] ;
				if ( !is_cutpoint_in_full_pose[ i_full ] && i_full < nres_full ){ // good, there's still a gap!

					// direct addition of single residue (or potentially, a block of residues starting with this residue)
					swa_moves.push_back( SWA_Move( i_full + 1, Attachment( i_full, BOND_TO_PREVIOUS ), ADD ) );

					// bulge skip...
					if ( allow_skip_bulge_ &&
							 (i == nres || res_list[ i ] + 2 < res_list[ i+1 ]) &&
							 i_full < (nres_full - 1)  &&
							 !is_cutpoint_in_full_pose[ i_full + 1 ] ){
						// for now can only handle single residue additions via skip-bulge ("floating base")
						Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( i_full + 2 );
						Size const other_pose_idx_intervening = full_model_info.get_idx_for_other_pose_with_residue( i_full + 1 );
						if ( other_pose_idx == 0 && other_pose_idx_intervening == 0 )
							swa_moves.push_back( SWA_Move( i_full + 2, Attachment( i_full, JUMP_TO_PREV_IN_CHAIN ), ADD ) );
					}
				}

			}
		}

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == 1 ) ||
					 ( fold_tree.is_cutpoint( i-1 ) && (res_list[ i ]-1 > res_list[ i-1 ]) ) ) { // could be a 5' chain terminus

				Size const i_full = res_list[ i ];
				if ( i_full > 1 && !is_cutpoint_in_full_pose[ i_full - 1 ] ) { // good, there's still a gap!

					// direct addition of single residue (or potentially, a block of residues starting with this residue)
					swa_moves.push_back( SWA_Move( i_full - 1, Attachment( i_full, BOND_TO_NEXT), ADD ) );

					// bulge skip...
					if ( allow_skip_bulge_ &&
							 (i == 1 || res_list[ i ] - 2 > res_list[ i - 1 ]) &&
							 i_full > 2 &&
							 !is_cutpoint_in_full_pose[ i_full - 2 ] ){
						// for now can only handle single residue additions via skip-bulge ("floating base")
						Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( i_full - 2 );
						Size const other_pose_idx_intervening = full_model_info.get_idx_for_other_pose_with_residue( i_full - 1 );
						if ( other_pose_idx == 0  && other_pose_idx_intervening == 0)
							swa_moves.push_back( SWA_Move( i_full - 2, Attachment( i_full, JUMP_TO_NEXT_IN_CHAIN ), ADD ) );
					}
				}

			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intramolecular_delete_move_elements( pose::Pose const & pose,
																														 utility::vector1< SWA_Move > & swa_moves ) {

		utility::vector1< SWA_Move > swa_moves_delete;
		get_intramolecular_split_move_elements( pose, swa_moves_delete, DELETE );

		// don't delete the only pose left.
		if ( const_full_model_info( pose ).other_pose_list().size() == 0 && pose.total_residue() == 2 ) return;

		for ( Size n = 1; n <= swa_moves_delete.size(); n++ ) swa_moves.push_back( swa_moves_delete[ n ] );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intramolecular_split_move_elements( pose::Pose const & pose,
																														utility::vector1< SWA_Move > & swa_moves,
																														MoveType const move_type ) {

		using namespace protocols::stepwise::sampling::rna;

		utility::vector1< SWA_Move > swa_moves_split;
		utility::vector1< bool > partition_definition;
		utility::vector1< Size > partition_res1, partition_res2;
		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const chains = figure_out_chains_from_full_model_info_const( pose );

		// first look at suites
		for ( Size i = 1; i < pose.total_residue(); i++ ){
			if ( !pose.fold_tree().is_cutpoint( i ) &&
					 ( domain_map[ i ] != domain_map[ i+1 ] ||
						 domain_map[ i ] == 0 || domain_map [ i+1 ] == 0 ) ){

				partition_definition = get_partition_definition( pose, i /*suite*/ );
				if ( partition_splits_a_fixed_domain( partition_definition, domain_map ) ) continue;

				partition_res1 = get_partition_res( partition_definition, true );
				partition_res2 = get_partition_res( partition_definition, false );

				if ( check_for_fixed_domain_or_from_scratch( pose, partition_res1 )  ) {
					swa_moves_split.push_back( SWA_Move( full_model_info.sub_to_full(partition_res2),
																							 Attachment( full_model_info.sub_to_full(i), BOND_TO_PREVIOUS ), move_type ) );
				}
				if ( check_for_fixed_domain_or_from_scratch( pose, partition_res2 ) ) {
					swa_moves_split.push_back( SWA_Move( full_model_info.sub_to_full(partition_res1),
																							 Attachment( full_model_info.sub_to_full(i+1), BOND_TO_NEXT ), move_type ) );
				}

			}
		}

		// now look at jumps. Note that, currently, additions by jump are only for single residues, so only
		// split cases in which single residues are deleted.
		for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){

			partition_definition = get_partition_definition_by_jump( pose, n );
			if ( partition_splits_a_fixed_domain( partition_definition, domain_map ) ) continue;

			Size const moving_res = pose.fold_tree().downstream_jump_residue( n );
			Size const reference_res = pose.fold_tree().upstream_jump_residue( n );
			if ( chains[ moving_res ] != chains[ reference_res ] ) continue;

			partition_res1 = get_partition_res( partition_definition, true );
			partition_res2 = get_partition_res( partition_definition, false );
			if ( partition_res1.size() == 1 && check_for_fixed_domain_or_from_scratch( pose, partition_res2 ) ){
				Size const anchor_res = get_anchor_res( partition_res1[1], pose );
				AttachmentType type = ( anchor_res > partition_res1[1] ) ? JUMP_TO_NEXT_IN_CHAIN : JUMP_TO_PREV_IN_CHAIN;
				swa_moves_split.push_back( SWA_Move( full_model_info.sub_to_full(partition_res1),
																						 Attachment( full_model_info.sub_to_full(anchor_res), type ), move_type ) );
			}
			if ( partition_res2.size() == 1 && check_for_fixed_domain_or_from_scratch( pose, partition_res1 ) ){
				Size const anchor_res = get_anchor_res( partition_res2[1], pose );
				AttachmentType type = ( anchor_res > partition_res2[1] ) ? JUMP_TO_NEXT_IN_CHAIN : JUMP_TO_PREV_IN_CHAIN;
				swa_moves_split.push_back( SWA_Move( full_model_info.sub_to_full(partition_res2),
																						 Attachment( full_model_info.sub_to_full(anchor_res), type ), move_type ) );
			}
		}
		for ( Size n = 1; n <= swa_moves_split.size(); n++ ) swa_moves.push_back( swa_moves_split[ n ] );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// intermolecular (cross-chain)
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intermolecular_split_move_elements( pose::Pose const & pose,
																														utility::vector1< SWA_Move > & swa_moves,
																														MoveType const move_type ){
		using namespace protocols::stepwise::sampling::rna;

		utility::vector1< SWA_Move > swa_moves_split;
		utility::vector1< bool > partition_definition;
		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const chains = figure_out_chains_from_full_model_info_const( pose );

		// look at jumps.
		for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			partition_definition = get_partition_definition_by_jump( pose, n );
			if ( partition_splits_a_fixed_domain( partition_definition, domain_map ) ) continue;

			Size const moving_res = pose.fold_tree().downstream_jump_residue( n );
			Size const reference_res = pose.fold_tree().upstream_jump_residue( n );
			if ( chains[ moving_res ] == chains[ reference_res ] ) continue;

			utility::vector1< Size > moving_partition_res = get_partition_res( partition_definition,  partition_definition[ moving_res ] );
			swa_moves_split.push_back( SWA_Move( full_model_info.sub_to_full( moving_partition_res ),
																					 Attachment( full_model_info.sub_to_full(reference_res), JUMP_INTERCHAIN ), move_type ) );
		}

		for ( Size n = 1; n <= swa_moves_split.size(); n++ ) swa_moves.push_back( swa_moves_split[ n ] );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	SWA_MoveSelector::check_for_fixed_domain_or_from_scratch(  pose::Pose const & pose,
																														 utility::vector1< Size> const & partition_res ) const {
		if ( resampling_ ) return true;
		if ( from_scratch_frequency_ > 0.0 ) return true;
		return check_for_fixed_domain( pose, partition_res );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// The point of this function is to allow big move_elements to remain fixed. It may not be entirely necessary.
	void
	SWA_MoveSelector::remove_from_consideration_first_multi_residue_move_element( utility::vector1< SWA_Move > & swa_moves,
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

		if ( remove_even_if_not_singlet || (number_of_multi_residue_move_elements == 1) ){

			utility::vector1< SWA_Move > swa_moves_new;

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
	SWA_MoveSelector::get_from_scratch_add_and_delete_elements( pose::Pose const & pose,
																															utility::vector1< SWA_Move > & swa_moves,
																															utility::vector1< Size > const & sample_res /*leave empty if no filter*/) {

		get_from_scratch_add_move_elements( pose, swa_moves );
		get_from_scratch_delete_move_elements( pose, swa_moves );
		if ( sample_res.size() > 0 ) filter_by_sample_res( swa_moves, sample_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_from_scratch_add_move_elements( pose::Pose const & pose,
																												utility::vector1< SWA_Move > & swa_moves ){
		// for RNA, these are dinucleotides that would be 'freely floating' -- not attached
		// to the current pose (or its sisters).
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		Size nres_full = full_model_info.full_sequence().size();

		Attachments const blank_attachments;

		for ( Size n = 1; n < nres_full; n++ ){
			if ( cutpoint_open_in_full_model.has_value( n ) ) continue; // must be contiguous dinucleotide.
			if ( already_instantiated_in_pose( pose, n     ) ) continue;
			if ( already_instantiated_in_pose( pose, n + 1 ) ) continue;
			// if ( n > 1 &&
			// 		 !cutpoint_open_in_full_model.has_value( n - 1 ) &&
			// 		 already_instantiated_in_pose( pose, n - 1 ) ) continue;
			// if ( (n + 1) < nres_full &&
			// 		 !cutpoint_open_in_full_model.has_value( n + 1 ) &&
			// 		 already_instantiated_in_pose( pose, n + 2 ) ) continue;
			swa_moves.push_back( SWA_Move( utility::tools::make_vector1( n, n+1 ), blank_attachments, FROM_SCRATCH ) );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::figure_out_from_scratch_delete_move_elements( pose::Pose const & pose,
																																	utility::vector1< SWA_Move > & swa_moves ) const {

		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		if ( pose.total_residue() == 2 &&
				 !partition_splits_a_fixed_domain( utility::tools::make_vector1(true,false), domain_map ) ) {
			swa_moves.push_back( SWA_Move( sub_to_full( 1, pose ), Attachment( sub_to_full( 2, pose ), BOND_TO_NEXT ), DELETE ) );
			swa_moves.push_back( SWA_Move( sub_to_full( 2, pose ), Attachment( sub_to_full( 1, pose ), BOND_TO_PREVIOUS ), DELETE ) );
		}
		utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
		for ( Size k = 1; k <= other_pose_list.size(); k++ ){
			figure_out_from_scratch_delete_move_elements( *( other_pose_list[k] ), swa_moves );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_from_scratch_delete_move_elements( pose::Pose const & pose,
																													 utility::vector1< SWA_Move > & swa_moves ){
		utility::vector1< SWA_Move > swa_moves_delete;
		figure_out_from_scratch_delete_move_elements( pose, swa_moves_delete );

		// don't delete the only pose left.
		if ( const_full_model_info( pose ).other_pose_list().size() == 0 && pose.total_residue() == 2 ) return;

		for ( Size n = 1; n <= swa_moves_delete.size(); n++ ) swa_moves.push_back( swa_moves_delete[n] );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intermolecular_add_and_delete_elements( pose::Pose const & pose,
																																utility::vector1< SWA_Move > & swa_moves,
																																utility::vector1< Size > const & sample_res /*leave empty if no filter*/) {

		get_intermolecular_add_move_elements( pose, swa_moves );
		get_intermolecular_delete_move_elements( pose, swa_moves );
		if ( sample_res.size() > 0 ) filter_by_sample_res( swa_moves, sample_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intermolecular_add_move_elements( pose::Pose const & pose,
																													utility::vector1< SWA_Move > & swa_moves ){

		using namespace core::pose::full_model_info;

		Size const & nres( pose.total_residue() );
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
		utility::vector1< Size > const chains_current = figure_out_chains_from_full_model_info_const( pose );
		utility::vector1< Size > const chains_full = get_chains_full( pose );
		utility::vector1< Size > const current_unique_chains = get_unique_chains( pose );

		if ( res_list.size() == 0 ) return;

		for ( Size i = 1; i <= nres; i++ ){
			Size const chain_i = chains_current[ i ];
			Size const i_full = res_list[ i ] ;

			if ( only_dock_preexisting_chunks_ ){
				for ( Size q = 1; q <= other_pose_list.size(); q++ ){
					Pose const & other_pose = *other_pose_list[ q ];
					for ( Size j = 1; j <= other_pose.total_residue(); j++ ){
						Size const j_full = get_res_list_from_full_model_info_const( other_pose )[ j ];
						if ( !allow_shared_chains_in_dock_poses_ && share_chains( current_unique_chains, pose, j_full ) ) break;
						Size const chain_j = chains_full[ j_full ];
						if ( chain_j == chain_i  ) continue;
						swa_moves.push_back( SWA_Move( j_full, Attachment( i_full, JUMP_INTERCHAIN ), ADD ) );
					}
				}
			} else {
				for ( Size j_full = 1; j_full <= full_model_info.size(); j_full++ ){
					if ( res_list.has_value( j_full ) ) continue;
					Size const chain_j = chains_full[ j_full ];
					if ( chain_j == chain_i  ) continue;
					if ( !allow_shared_chains_in_dock_poses_ && share_chains( current_unique_chains, pose, j_full ) ) continue;
					swa_moves.push_back( SWA_Move( j_full, Attachment( i_full, JUMP_INTERCHAIN ), ADD ) );
				}

			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	SWA_MoveSelector::get_unique_chains( pose::Pose const & pose ){
		utility::vector1< Size > const chains = figure_out_chains_from_full_model_info_const( pose );
		utility::vector1< Size > unique_chains;
		for ( Size n = 1; n <= chains.size(); n++ ){
			if ( !unique_chains.has_value( chains[n] ) ) unique_chains.push_back( chains[n] );
		}
		return unique_chains;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// this is pretty slow, actually -- can be sped up by pre-computing chains_full, unique chains in each pose, etc.
	bool
	SWA_MoveSelector::share_chains( utility::vector1< Size > const & current_unique_chains,
																	pose::Pose const & pose,
																	Size const j_full ) {

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const chains_full = get_chains_full( pose );

		Size const chain_j = chains_full[ j_full ];
		utility::vector1< Size > other_pose_chains = utility::tools::make_vector1( chain_j );

		Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( j_full );
		if ( other_pose_idx > 0 ) other_pose_chains = get_unique_chains( *full_model_info.other_pose_list()[ other_pose_idx ] );

		for ( Size q = 1; q <= current_unique_chains.size(); q++ ){
			if ( other_pose_chains.has_value( current_unique_chains[ q ] ) ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_intermolecular_delete_move_elements( pose::Pose const & pose,
																														 utility::vector1< SWA_Move > & swa_moves ){
		get_intermolecular_split_move_elements( pose, swa_moves, DELETE );
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::filter_by_sample_res( utility::vector1< SWA_Move > & swa_moves,
																					utility::vector1< Size > const & sample_res ){

		if ( sample_res.size() == 0 ) return;

		utility::vector1< SWA_Move > swa_moves_filtered;

		for ( Size i = 1; i <= swa_moves.size(); i++ ){

			SWA_Move const & swa_move = swa_moves[ i ];

			if ( swa_move.move_type() == ADD || swa_move.move_type() == FROM_SCRATCH ){
				MoveElement add_element = swa_move.move_element();
				bool in_sample_res( true );
				for ( Size n = 1; n <= add_element.size(); n++ ){
					if ( !sample_res.has_value( add_element[ n ] ) ) {
						in_sample_res = false;
						break;
					}
				}
				if ( !in_sample_res ) continue;
			}
			swa_moves_filtered.push_back( swa_move );
		}

		swa_moves = swa_moves_filtered;
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// helper functions
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// don't break up any fixed domains by a split.
	bool
	SWA_MoveSelector::partition_splits_a_fixed_domain( utility::vector1< Size > const & partition_definition,
																										 utility::vector1< Size > const & domain_map ) const {
		// parse out which residues go into each domain.
		utility::vector1< Size > blank_vector;
		utility::vector1< utility::vector1< Size > > all_domain_res;
		for ( Size i = 1; i <= domain_map.size(); i++ ){
			Size const & domain = domain_map[ i ];
			for ( Size n = all_domain_res.size(); n < domain; n++ ) all_domain_res.push_back( blank_vector );
			if ( domain > 0 ) all_domain_res[ domain ].push_back( i );
		}
		for ( Size n = 1; n <= all_domain_res.size(); n++ ){
			// for this domain, check if all residues are in a single partition.
			bool in_partition_0( false ), in_partition_1( false );
			for ( Size k = 1; k <= all_domain_res[ n ].size(); k++ ){
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
	SWA_MoveSelector::get_terminal_move_elements( pose::Pose const & pose,
																								utility::vector1< SWA_Move > & swa_moves,
																								MoveType const & move_type ) {

		using namespace core::pose::full_model_info;

		if ( get_res_list_from_full_model_info_const( pose ).size() == 0 ) return;
		utility::vector1< MoveElement > const & move_elements = get_move_elements_from_full_model_info_const( pose );

		Size const num_move_elements( move_elements.size() );
		for ( Size n = 1; n <= num_move_elements; n++ ){

			MoveElement const & move_element = move_elements[ n ];
			utility::vector1< Attachment > attachments = get_attachments( pose, move_element );
			// at least for now, each move_element should be connected to another one.
			if ( attachments.size() == 0 ) {
				runtime_assert( num_move_elements == 1 );
				continue;
			}
			if ( attachments.size() > 1 ) continue; // not a terminal

			// for now, floating base can only handle single residues.
			if ( attachments.size() == 1 &&
					 ( attachments[1].attachment_type() == JUMP_TO_PREV_IN_CHAIN || attachments[1].attachment_type() == JUMP_TO_NEXT_IN_CHAIN ) &&
					 move_element.size() > 1 ) continue;

			swa_moves.push_back( SWA_Move( move_element, attachments, move_type ) );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	SWA_MoveSelector::get_internal_move_elements( pose::Pose const & pose,
																								utility::vector1< SWA_Move > & swa_moves,
																								MoveType const & move_type ) {

		using namespace core::pose::full_model_info;

		utility::vector1< MoveElement > const & move_elements = get_move_elements_from_full_model_info_const( pose );

		Size const num_move_elements( move_elements.size() );
		for ( Size n = 1; n <= num_move_elements; n++ ){

			MoveElement const & move_element = move_elements[ n ];
			utility::vector1< Attachment > attachments = get_attachments( pose, move_element );

			// at least for now, each move_element should be connected to another one.
			if ( attachments.size() != 2 ) continue;

			// note that PREVIOUS should be before NEXT.
			if ( attachments[1].attachment_type() == BOND_TO_PREVIOUS &&
					 attachments[2].attachment_type() == BOND_TO_NEXT ) {
				swa_moves.push_back( SWA_Move( move_element, attachments, move_type ) );
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	SWA_MoveSelector::get_partition_res( utility::vector1< bool > const & partition_definition, bool const setting ) const {
		utility::vector1< Size > partition_res;
		for ( Size n = 1; n <= partition_definition.size(); n++ ){
			if ( partition_definition[ n ]  == setting ) partition_res.push_back( n );
		}
		return partition_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Attachments
	SWA_MoveSelector::get_attachments( pose::Pose const & pose, Size const & moving_res ){
		return get_attachments( pose, utility::tools::make_vector1( moving_res ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// A moving element (or even a single residue ) can actually have several connections.
	// Some are covalent (polymeric) connections (BOND_TO_PREVIOUS, BOND_TO_NEXT)
	// Some could involve jumps to residues in same chain (e.g., in 'skip-bulge' moves).
	// Some could involve jumps to other chains or ligands.
	Attachments
	SWA_MoveSelector::get_attachments( pose::Pose const & pose, MoveElement const & move_element ){

		using namespace core::pose::full_model_info;

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );

		Attachments attachments_to_previous, attachments_to_next, attachments;

		// need to find number of attachment points of this move_element to rest of pose.
		for ( Size k = 1; k <= move_element.size(); k++ ){

			Size const i_full = move_element[ k ];
			Size const i = res_list.index( i_full );

			if ( i > 1 && ! move_element.has_value( res_list[ i - 1 ] ) ){
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
				if ( !fold_tree.is_cutpoint( i ) ){
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
		//		runtime_assert( attachments.size() <= 4 );
		return attachments;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	SWA_MoveSelector::already_instantiated_in_pose( pose::Pose const & pose, Size const & resnum_in_full_model_numbering ) const {

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		if ( res_list.size() == 0 ) return false;

		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( res_list[ n ] == resnum_in_full_model_numbering ) return true;
		}

		utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
		for ( Size k = 1; k <= other_pose_list.size(); k++ ){
			if ( already_instantiated_in_pose( *(other_pose_list[ k ]), resnum_in_full_model_numbering ) ) return true;
		}

		return false;
	}


} //monte_carlo
} //stepwise
} //protocols
