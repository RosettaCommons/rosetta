// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/CartesianMultifunc.hh
/// @brief
/// @author Frank DiMaio

/// Unit headers
#include <core/optimization/CartesianMultifunc.hh>

/// Package headers
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/cartesian_minimize.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

/// Project headers
#include <basic/prof.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Energies.hh>

/// Utility headers
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {

Real
CartesianMultifunc::operator ()( Multivec const & vars ) const {
	PROF_START( basic::FUNC );
	min_map_.copy_dofs_to_pose( pose_, vars );
	Real const score( score_function_( pose_ ) );
	PROF_STOP( basic::FUNC );
	return score;
}

void
CartesianMultifunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const
{
	PROF_START( basic::DFUNC );
	// in cartesian_minimize.cc
	cartesian_dfunc( pose_, min_map_, score_function_, vars, dE_dvars );
	// optional derivative checking
	if ( deriv_check_ ) {
		cart_numerical_derivative_check( min_map_, *this, vars, dE_dvars, deriv_check_result_, deriv_check_verbose_ );
	}
	PROF_STOP( basic::DFUNC );
}

void CartesianMultifunc::set_deriv_check_result( NumericalDerivCheckResultOP deriv_check_result )
{
	deriv_check_result_ = deriv_check_result;
}

core::pose::Pose & CartesianMultifunc::pose() const {
	return pose_;
}

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this function.
void
CartesianMultifunc::dump( Multivec const & vars, Multivec const & vars2 ) const {
    // AMW: cppcheck flags these bools as being inappropriately scoped for obvious reasons
	// I don't suggest changing them, because the point is to have them at the top	
	bool debug_inaccurateG = false; // disables everything below
	bool check_score_components = true;
	bool check_score_components_verbose = false;
	bool check_rama = false;
	bool check_hbonds = true;
	//bool check_nblist = true;

	if ( ! debug_inaccurateG ) return;

	static int count_dumped( 0 );
	static bool after( true ); // dump two poses, a before and an after. Note, dumping poses as pdbs is often useless.

	if ( after ) { ++count_dumped; after = false; }
	else { after = true; }

	pose::Pose pose1( pose_ );
	pose::Pose pose2( pose_ );
	min_map_.copy_dofs_to_pose( pose1, vars );
	min_map_.copy_dofs_to_pose( pose2, vars2 );

	min_map_.copy_dofs_to_pose( pose_, vars );
	Real score_vars( score_function_( pose_ ) );

	min_map_.copy_dofs_to_pose( pose_, vars2 );
	Real score_vars2( score_function_( pose_ ) );

	Real alt_score_vars = score_function_( pose1 );
	pose1.dump_pdb( "Cartesian_multifunc_error_pose_before" + utility::to_string( count_dumped  ) + ".pdb" );

	Real alt_score_vars2 = score_function_( pose2 );
	pose2.dump_pdb( "Cartesian_multifunc_error_pose_after" + utility::to_string( count_dumped  ) + ".pdb" );

	std::cerr << "starting pose energies: " << score_vars << " vs " << alt_score_vars << std::endl;
	pose1.energies().total_energies().show_weighted( std::cerr, score_function_.weights() );
	std::cerr << std::endl;
	std::cerr << "moved pose energies: " << score_vars2 << " vs " << alt_score_vars2 << std::endl;
	pose2.energies().total_energies().show_weighted( std::cerr, score_function_.weights() );
	std::cerr << std::endl;
	using namespace scoring;

	if ( check_score_components ) {
		// slow! Iterate through all the components and check their derivatives one by one.
		const_cast< bool & > (deriv_check_) = true;
		if ( check_score_components_verbose ) {
			const_cast< bool & > (deriv_check_verbose_) = true;
		}
		Multivec dvars( vars );
		scoring::EnergyMap orig_weights( score_function_.weights() );
		for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
			using namespace scoring;

			if ( score_function_.weights()[ (ScoreType ) ii ] == 0.0 ) continue;

			for ( Size jj = 1; jj <= scoring::n_score_types; ++jj ) {
				if ( jj == ii ) {
					const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType) jj, orig_weights[ (scoring::ScoreType) jj ]);
				} else if ( score_function_.weights()[ (scoring::ScoreType ) jj ] != 0.0 ) {
					const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType) jj, 1e-9 );
				}
			}
			std::cout << "Checking score type: " << scoring::ScoreType( ii ) << std::endl;
			dfunc( vars, dvars ); // invokes numeric derivative checker.
		}
		for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
			if ( orig_weights[ scoring::ScoreType( ii ) ] != 0 ) {
				const_cast< scoring::ScoreFunction & > (score_function_).set_weight( (scoring::ScoreType)ii, orig_weights[ (scoring::ScoreType) ii ]);
			}
		}
		const_cast< bool & > (deriv_check_) = false;
		const_cast< bool & > (deriv_check_verbose_) = false;
	}

	if ( check_rama ) {
		// useful if rama seems to be the culprit.  This is the only piece of code
		// that invokes eval_rama_score_all, so that function may be hacked to provide
		// clearer debugging output.
		ScoringManager::get_instance()->get_Ramachandran().eval_rama_score_all( pose1, score_function_ );
	}

	if ( check_hbonds ) {
		scoring::hbonds::HBondSet hbond_set;
		fill_hbond_set( pose1, true, hbond_set );

		for ( Size ii = 1; ii <= hbond_set.nhbonds(); ++ii ) {
			scoring::hbonds::HBond const & iihbond = hbond_set.hbond( ii );
			std::cerr << "Hbond " << ii <<
				" d: " << iihbond.don_res() << " " << iihbond.don_hatm() <<
				" a: " << iihbond.acc_res() << " " << iihbond.acc_atm() <<
				" e: " << iihbond.energy() << " " << iihbond.weight() <<
				" good? " << hbond_set.allow_hbond( iihbond ) << std::endl;
		}
	}
}

CartesianMultifunc::CartesianMultifunc(
	pose::Pose & pose_in,
	CartesianMinimizerMap & min_map_in,
	scoring::ScoreFunction const & scorefxn_in,
	bool const deriv_check_in,
	bool const deriv_check_verbose_in
) :
	pose_( pose_in ),
	min_map_( min_map_in ),
	score_function_( scorefxn_in ),
	deriv_check_( deriv_check_in ),
	deriv_check_verbose_( deriv_check_verbose_in ),
	deriv_check_result_( /* 0 */ )
{}

CartesianMultifunc::~CartesianMultifunc() {}

CartesianMinimizerMap const & 
CartesianMultifunc::min_map() const {
	 return min_map_;
}

core::scoring::ScoreFunction const & CartesianMultifunc::score_function() const {
	 return score_function_;
}

} // namespace optimization
} // namespace core

