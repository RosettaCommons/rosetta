// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover
/// @brief Torsions an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/SWA_MonteCarloUtil.hh>
#include <protocols/swa/monte_carlo/SWA_Move.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/util.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

#include <map>

#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG(239111);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.swa.monte_carlo.SWA_MonteCarloUtil" ) ;

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// I think this is deprecated in favor of get_moving_chunk_case...
	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, Size const i ) {

		// is this in use?
		runtime_assert( 1 == 0 );

		MovingResidueCase moving_residue_case( NO_CASE );

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );
		if ( i == nres || fold_tree.is_cutpoint( i ) ){ // could be a 5' chain terminus
			if ( i == 1 || fold_tree.is_cutpoint( i-1 ) ){
				moving_residue_case = FLOATING_BASE; // don't know how to handle this yet.
			} else {
				moving_residue_case = CHAIN_TERMINUS_3PRIME;
			}
		} else if ( i == 1 || fold_tree.is_cutpoint( i-1 ) ){
			moving_residue_case = CHAIN_TERMINUS_5PRIME;
		} else {
			moving_residue_case = INTERNAL;
		}

		return moving_residue_case;
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_delete_chunks( pose::Pose & pose,
															 utility::vector1< SWA_Move > & swa_moves ) {

		swa_moves.clear();
		get_potential_terminal_chunks( pose, swa_moves, DELETE );

		// don't delete a multi_residue_chunk if its the only one!
		remove_from_consideration_first_multi_residue_chunk( swa_moves, false /*remove_even_if_not_singlet*/ );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// is following in use anymore?
	void
	get_potential_resample_chunks( pose::Pose & pose,
																 utility::vector1< Chunk > & possible_chunks ){

		utility::vector1< SWA_Move > swa_moves;
		get_potential_resample_chunks( pose, swa_moves );
		for ( Size n = 1; n <= swa_moves.size(); n++ ) possible_chunks.push_back( swa_moves[ n ].chunk() );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_resample_chunks( pose::Pose & pose,
																 utility::vector1< SWA_Move > & swa_moves ) {

		swa_moves.clear();
		get_potential_terminal_chunks( pose, swa_moves, NO_ADD_OR_DELETE );
		remove_from_consideration_first_multi_residue_chunk( swa_moves, true /*remove_even_if_not_singlet*/ );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// The point of this function is to allow big chunks to remain fixed. It may not be entirely necessary.
	void
	remove_from_consideration_first_multi_residue_chunk( utility::vector1< SWA_Move > & swa_moves,
																											 bool remove_even_if_not_singlet ){

		Size number_of_multi_residue_chunks( 0 ), first_multi_residue_chunk( 0 );
		for ( Size n = 1; n <= swa_moves.size(); n++ ) {
			if ( swa_moves[ n ].chunk().size() > 1 ) {
				number_of_multi_residue_chunks++;
				if ( first_multi_residue_chunk == 0 ) first_multi_residue_chunk = n;
			}
		}
		if ( number_of_multi_residue_chunks == 0 ) return;

		runtime_assert( number_of_multi_residue_chunks > 0 );

		if ( remove_even_if_not_singlet || (number_of_multi_residue_chunks == 1) ){

			utility::vector1< SWA_Move > swa_moves_new;

			for ( Size n = 1; n <= swa_moves.size(); n++ ) {
				if ( n == first_multi_residue_chunk ) continue;
				swa_moves_new.push_back( swa_moves[ n ] );
			}
			swa_moves = swa_moves_new;

		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_terminal_chunks( pose::Pose & pose,
																 utility::vector1< SWA_Move > & swa_moves,
																 AddOrDeleteChoice const & choice ) {

		using namespace core::pose::full_model_info;

		swa_moves.clear();

		utility::vector1< Chunk > const & moving_chunks = get_moving_chunks_from_full_model_info( pose );

		Size const num_chunks( moving_chunks.size() );

		for ( Size n = 1; n <= num_chunks; n++ ){

			Chunk const & moving_chunk = moving_chunks[ n ];
			MovingResidueCase const moving_residue_case = get_moving_residue_case( pose, moving_chunk );

			if ( moving_residue_case == INTERNAL ) continue;
			if ( moving_residue_case == NO_CASE ) runtime_assert( num_chunks == 1 );

			swa_moves.push_back( SWA_Move( moving_chunk, moving_residue_case, choice ) );
		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, utility::vector1< Size > const & moving_chunk ){

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );

		// need to find number of attachment points of this chunk to rest of pose.
		bool attached_on_five_prime( false ), attached_on_three_prime( false );
		for ( Size k = 1; k <= moving_chunk.size(); k++ ){

			Size const i = moving_chunk[ k ];

			// if at the ends of the pose, certainly no attachment.
			if ( i == 1 ) continue;
			if ( i == nres ) continue;

			// potential connection to previous residue in rest of pose
			if ( ! moving_chunk.has_value( i - 1 ) ){
				if ( !fold_tree.is_cutpoint( i - 1 ) )  {
					runtime_assert( !attached_on_five_prime ); // can't have more than one five_prime attachment ?
					attached_on_five_prime = true;
				}
			}

			// potential connection to next residue in rest of pose
			if ( ! moving_chunk.has_value( i + 1 ) ){
				runtime_assert( !attached_on_three_prime );  // can't have more than one three_prime attachment ?
				if ( !fold_tree.is_cutpoint( i ) )  attached_on_three_prime = true;
			}
		}

		if ( attached_on_five_prime && attached_on_three_prime ){
			return INTERNAL;
		} else if ( !attached_on_five_prime && attached_on_three_prime ){
			return CHAIN_TERMINUS_5PRIME;
		} else if ( !attached_on_three_prime && attached_on_five_prime ){
			return CHAIN_TERMINUS_3PRIME;
		}
		return NO_CASE;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_add_chunks( pose::Pose & pose,
														utility::vector1< SWA_Move > & swa_moves ){

		using namespace core::pose::full_model_info;

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );

		FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		Size nres_full = full_model_info.full_sequence().size();

		utility::vector1< bool > is_cutpoint_in_full_pose;
		for ( Size i = 1; i <= nres_full; i++ ) is_cutpoint_in_full_pose.push_back( false );
		for ( Size n = 1; n <= cutpoint_open_in_full_model.size(); n++ ) is_cutpoint_in_full_pose[ cutpoint_open_in_full_model[n] ] = true;

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == nres ) ||
					 ( fold_tree.is_cutpoint( i ) && (res_list[ i ]+1 < res_list[ i+1 ]) ) ) { // could be a 3' chain terminus

				Size const i_full = res_list[ i ] ;
				if ( !is_cutpoint_in_full_pose[ i_full ] && i_full < nres_full ){ // good, there's still a gap!
					swa_moves.push_back( SWA_Move( utility::tools::make_vector1( i ), CHAIN_TERMINUS_3PRIME, ADD ) );
				}

			}
		}

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == 1 ) ||
					 ( fold_tree.is_cutpoint( i-1 ) && (res_list[ i ]-1 > res_list[ i-1 ]) ) ) { // could be a 5' chain terminus

				Size const i_full = res_list[ i ];
				if ( i_full > 1 && !is_cutpoint_in_full_pose[ i_full-1 ] ) { // good, there's still a gap!
					swa_moves.push_back( SWA_Move( utility::tools::make_vector1( i ), CHAIN_TERMINUS_5PRIME, ADD ) );
				}

			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_random_chunk_at_chain_terminus( pose::Pose & pose,
																			SWA_Move & swa_move,
																			bool const disallow_delete,
																			bool const disallow_resample ){

		utility::vector1< Size > sample_res; /*leave empty if no filter*/
		get_random_chunk_at_chain_terminus( pose, swa_move,
																				disallow_delete, disallow_resample,
																				sample_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_random_chunk_at_chain_terminus( pose::Pose & pose,
																			SWA_Move & swa_move,
																			bool const disallow_delete,
																			bool const disallow_resample,
																			utility::vector1< Size > const & sample_res /*leave empty if no filter*/) {

		using namespace core::pose::full_model_info;

		utility::vector1< SWA_Move >  swa_moves;

		if ( !disallow_resample )  get_potential_resample_chunks( pose, swa_moves );
		if ( !disallow_delete )    get_potential_delete_chunks( pose, swa_moves );

		get_potential_add_chunks( pose, swa_moves );

		if ( sample_res.size() > 0 ) filter_by_sample_res( swa_moves,
																											 sample_res, get_res_list_from_full_model_info( pose ) );
		for ( Size n = 1; n <= swa_moves.size(); n++ )	TR.Debug << swa_moves[ n ] << std::endl;

		if ( swa_moves.size() == 0 ){
			utility::vector1 < Size > blank; //blank
			swa_move = SWA_Move( blank, NO_CASE, NO_ADD_OR_DELETE );
			return;
		}

		swa_move =  RG.random_element( swa_moves );

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	filter_by_sample_res( utility::vector1< SWA_Move > & swa_moves,
											  utility::vector1< Size > const & sample_res,
											  utility::vector1< Size > const & res_list ){

		if ( sample_res.size() == 0 ) return;

		utility::vector1< SWA_Move > swa_moves_filtered;

		for ( Size i = 1; i <= swa_moves.size(); i++ ){

			SWA_Move const & swa_move = swa_moves[ i ];

			if ( swa_move.add_or_delete_choice() == ADD ){
				runtime_assert( swa_move.chunk().size() == 1 );
				Size const & possible_res = swa_move.chunk() [ 1 ];

				// need to check if added residue is permitted, given sample res.
				if ( swa_move.moving_residue_case()  == CHAIN_TERMINUS_3PRIME ) { // prepend
					if ( ! sample_res.has_value( res_list[ possible_res ] + 1 ) ) continue;
				} else {
					runtime_assert( swa_move.moving_residue_case()  == CHAIN_TERMINUS_5PRIME ) // prepend
					if ( ! sample_res.has_value( res_list[ possible_res ] - 1 ) ) continue;
				}
			}

			swa_moves_filtered.push_back( swa_move );
		}

		swa_moves = swa_moves_filtered;
	}


}
}
}
