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

#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>
#include <protocols/swa/monte_carlo/types.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/moves/MonteCarlo.hh>

#include <basic/Tracer.hh>

#include <map>

#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG(239111);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_swa_monte_carlo_util" ) ;

using namespace core;
using core::Real;

namespace protocols {
namespace swa {
namespace monte_carlo {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MovingResidueCase
	get_moving_residue_case( pose::Pose const & pose, Size const i ) {

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
	get_potential_delete_residues( pose::Pose & pose,
																 utility::vector1< Size > & possible_res,
																 utility::vector1< MovingResidueCase > & moving_residue_cases,
																 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices ) {


		get_potential_terminal_residues( pose, possible_res, moving_residue_cases, add_or_delete_choices, DELETE );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_resample_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res ){

		utility::vector1< MovingResidueCase > moving_residue_cases;
		utility::vector1< AddOrDeleteChoice > add_or_delete_choices;
		get_potential_terminal_residues( pose, possible_res, moving_residue_cases, add_or_delete_choices, NO_ADD_OR_DELETE );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_resample_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res,
																	 utility::vector1< MovingResidueCase > & moving_residue_cases,
																	 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices ) {

		get_potential_terminal_residues( pose, possible_res, moving_residue_cases, add_or_delete_choices, NO_ADD_OR_DELETE );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_terminal_residues( pose::Pose & pose,
																	 utility::vector1< Size > & possible_res,
																	 utility::vector1< MovingResidueCase > & moving_residue_cases,
																	 utility::vector1< AddOrDeleteChoice > & add_or_delete_choices,
																	 AddOrDeleteChoice const & choice ) {


		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );

		SubToFullInfo & sub_to_full_info = nonconst_sub_to_full_info_from_pose( pose );
		utility::vector1< Size > const & moving_res_list = sub_to_full_info.moving_res_list();

		for ( Size n = 1; n <= moving_res_list.size(); n++ ){

			Size const i = moving_res_list[ n ];

			if ( i == nres || fold_tree.is_cutpoint( i ) ){ // could be a 3' chain terminus

				possible_res.push_back( i );
				moving_residue_cases.push_back( CHAIN_TERMINUS_3PRIME );
				add_or_delete_choices.push_back( choice );

			} else if ( i == 1 || fold_tree.is_cutpoint( i-1 ) ) {

				possible_res.push_back( i );
				moving_residue_cases.push_back( CHAIN_TERMINUS_5PRIME );
				add_or_delete_choices.push_back( choice );

			}


		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_potential_add_residues( pose::Pose & pose,
															utility::vector1< Size > & possible_res,
															utility::vector1< MovingResidueCase > & moving_residue_cases,
															utility::vector1< AddOrDeleteChoice > & add_or_delete_choices ) {

		Size const & nres( pose.total_residue() );
		kinematics::FoldTree const & fold_tree( pose.fold_tree() );

		SubToFullInfo & sub_to_full_info = nonconst_sub_to_full_info_from_pose( pose );
		std::map< Size, Size > sub_to_full = sub_to_full_info.sub_to_full();
		utility::vector1< Size > cutpoints_in_full_pose = sub_to_full_info.cutpoints_in_full_pose();
		Size nres_full = sub_to_full_info.full_sequence().size();

		utility::vector1< bool > is_cutpoint_in_full_pose;
		for ( Size i = 1; i <= nres_full; i++ ) is_cutpoint_in_full_pose.push_back( false );
		for ( Size n = 1; n <= cutpoints_in_full_pose.size(); n++ ) is_cutpoint_in_full_pose[ cutpoints_in_full_pose[n] ] = true;

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == nres ) ||
					 ( fold_tree.is_cutpoint( i ) && (sub_to_full[ i ]+1 < sub_to_full[ i+1 ]) ) ) { // could be a 3' chain terminus

				Size const i_full = sub_to_full[ i ] ;
				if ( !is_cutpoint_in_full_pose[ i_full ] && i_full < nres_full ){ // good, there's still a gap!

					possible_res.push_back( i );
					moving_residue_cases.push_back( CHAIN_TERMINUS_3PRIME );
					add_or_delete_choices.push_back( ADD );

				}
			}
		}

		for ( Size i = 1; i <= nres; i++ ){

			if ( ( i == 1 ) ||
					 ( fold_tree.is_cutpoint( i-1 ) && (sub_to_full[ i ]-1 > sub_to_full[ i-1 ]) ) ) { // could be a 5' chain terminus

				Size const i_full = sub_to_full[ i ];
				if ( i_full > 1 && !is_cutpoint_in_full_pose[ i_full-1 ] ) { // good, there's still a gap!

					possible_res.push_back( i );
					moving_residue_cases.push_back( CHAIN_TERMINUS_5PRIME );
					add_or_delete_choices.push_back( ADD );
				}
			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_random_residue_at_chain_terminus( pose::Pose & pose,
																				Size & residue_at_chain_terminus,
																				MovingResidueCase & moving_residue_case,
																				AddOrDeleteChoice & add_or_delete_choice,
																				bool const disallow_delete,
																				bool const disallow_resample ) {


		utility::vector1< Size >  possible_res;
		utility::vector1< MovingResidueCase > moving_residue_cases;
		utility::vector1< AddOrDeleteChoice > add_or_delete_choices;

		if ( !disallow_resample )  get_potential_resample_residues( pose,  possible_res, moving_residue_cases, add_or_delete_choices );

		if ( !disallow_delete )    get_potential_delete_residues( pose,  possible_res, moving_residue_cases, add_or_delete_choices );

		get_potential_add_residues( pose, possible_res, moving_residue_cases, add_or_delete_choices );

		Size const res_idx =  RG.random_range( 1, possible_res.size() );

		residue_at_chain_terminus = possible_res[ res_idx ];
		moving_residue_case       = moving_residue_cases[ res_idx ];
		add_or_delete_choice      = add_or_delete_choices[ res_idx ];

	}


}
}
}
