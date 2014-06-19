// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/checker/RNA_ChainClosureChecker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/sampling/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/id/TorsionID.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>

using namespace core;

static basic::Tracer TR( "protocols.stepwise.sampling.rna.checker.RNA_ChainClosureChecker" ) ;

///////////////////////////////////////////////////////////////////////////////////////////////////
// Kind of weird. I think chain closure always 'passes', since it depends on angle and atom-pair
//  constraints, which are not necessarily set up ahead of time.
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace checker {

	//Constructor
	RNA_ChainClosureChecker::RNA_ChainClosureChecker( pose::Pose const & pose, Size const five_prime_res ):
		chain_break_screening_pose_( pose ),
		five_prime_res_( five_prime_res ),
		reinitialize_CCD_torsions_( false ),
		verbose_( false )
	{
		chain_break_scorefxn_ =  new core::scoring::ScoreFunction;
		chain_break_scorefxn_->set_weight( core::scoring::angle_constraint, 1.0 );
		chain_break_scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	}

	//Destructor
	RNA_ChainClosureChecker::~RNA_ChainClosureChecker()
	{}


	////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChainClosureChecker::copy_CCD_torsions( pose::Pose & pose ) const {

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;

		Size const three_prime_res = five_prime_res_ + 1;
		//Even through there is the chain_break, alpha of 3' and epl and gamma of 5' should be defined due to the existence of the upper and lower variant type atoms.
 		copy_CCD_torsions_general( pose, five_prime_res_, three_prime_res );

	}


	////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChainClosureChecker::copy_CCD_torsions_general( pose::Pose & pose, Size const five_prime_res, Size const three_prime_res ) const {

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;

		if ( ( five_prime_res ) != ( three_prime_res - 1 ) ) utility_exit_with_message( "( five_prime_res ) != ( three_prime_res - 1 )" );

		conformation::Residue const & lower_res = chain_break_screening_pose_.residue( five_prime_res );
		conformation::Residue const & upper_res = chain_break_screening_pose_.residue( three_prime_res );

		for ( Size n = 1; n <= 3; n++ ){ //alpha, beta, gamma of 3' res
			pose.set_torsion( TorsionID( three_prime_res, id::BB,  n ), upper_res.mainchain_torsion( n ) );
		}

		for ( Size n = 5; n <= 6; n++ ){ //epsilon and zeta of 5' res
			pose.set_torsion( TorsionID( five_prime_res, id::BB,  n ), lower_res.mainchain_torsion( n ) );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChainClosureChecker::check_loop_closed( pose::Pose const & pose ){
		static protocols::farna::RNA_LoopCloser rna_loop_closer;
		return ( rna_loop_closer.check_closure( pose, five_prime_res_ ) );
	}

	////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChainClosureChecker::chain_break_screening_general( pose::Pose & chain_break_screening_pose,
																										 core::scoring::ScoreFunctionOP const & chain_break_scorefxn,
																										 Size const five_prime_res ){

		using namespace core::scoring;

		static protocols::farna::RNA_LoopCloser rna_loop_closer;
		runtime_assert( chain_break_screening_pose.residue( five_prime_res ).has_variant_type( chemical::CUTPOINT_LOWER ) );
		runtime_assert( chain_break_screening_pose.residue( five_prime_res + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) );

		if ( reinitialize_CCD_torsions_ ) set_CCD_torsions_to_zero( chain_break_screening_pose, five_prime_res );

		//		Real const mean_dist_err=rna_loop_closer.apply( chain_break_screening_pose, five_prime_res);
		rna_loop_closer.apply( chain_break_screening_pose, five_prime_res );

		( *chain_break_scorefxn )( chain_break_screening_pose );

		scoring::EMapVector & energy_map = chain_break_screening_pose.energies().total_energies();
		Real const angle_score = energy_map[scoring::angle_constraint];
		Real const distance_score = energy_map[scoring::atom_pair_constraint];

		if ( angle_score < 5 ) count_data_.good_angle_count++;
		if ( distance_score < 5 ) count_data_.good_distance_count++;
		if ( ( angle_score < 5 ) && ( distance_score < 5 ) ){
			count_data_.chain_break_screening_count++;
			if ( verbose_ ){
				//				TR.Debug << " C5_O3= " << C5_O3_distance << " C5_O3_n= " << count_data_.C5_O3_distance_count;
				TR.Debug << "  chain_closable_geometry_count = " << count_data_.chain_closable_geometry_count;
				TR.Debug << " angle = " << angle_score << " dist = " << distance_score;
				TR.Debug << " angle_n = " << count_data_.good_angle_count;
				TR.Debug << " dist_n = " << count_data_.good_distance_count;
				TR.Debug << " chain_break_screening = " << count_data_.chain_break_screening_count;
				TR.Debug << " tot = " << count_data_.tot_rotamer_count << std::endl;
			}
			return true;
		}

		return false;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChainClosureChecker::check_screen(){
		return chain_break_screening_general( chain_break_screening_pose_, chain_break_scorefxn_, five_prime_res_ );
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChainClosureChecker::check_screen( pose::Pose & pose ){
		return chain_break_screening_general( pose, chain_break_scorefxn_, five_prime_res_ );
	}


} //checker
} //rna
} //sampling
} //stepwise
} //protocols
