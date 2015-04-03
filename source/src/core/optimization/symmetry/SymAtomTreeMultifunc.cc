// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/SymAtomTreeMultifunc.hh
/// @brief
/// @author

/// Unit headers
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

/// Package headers
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/sym_atom_tree_minimize.hh>

/// Project headers


#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

/// Utility headers
#include <basic/prof.hh>

/// Numeric headers

#include <utility/vector1.hh>

//Auto Headers
#include <utility/string_util.hh>


namespace core {
namespace optimization {
namespace symmetry {

bool debug_inaccurateG = false;
bool check_score_components = false;
bool check_score_components_verbose = false;
bool check_rama = false;
bool check_hbonds = false;

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this file.
void
SymAtomTreeMultifunc::dump( Multivec const & vars, Multivec const & /*vars2*/ ) const
{
	if ( ! debug_inaccurateG ) return;

	static int count_dumped( 0 );
	static bool after( true ); // dump two poses, a before and an after. Note, dumping poses as pdbs is often useless.

	if ( after ) { ++count_dumped; after = false; }
	else { after = true; }

	symm_min_map_.copy_dofs_to_pose( pose_, vars );

	std::string which =  ( after ? "after_" : "before_" );
	pose_.dump_pdb( "atomtree_multifunc_error_pose_" + which + utility::to_string( count_dumped  ) + ".pdb" );
	score_function_( pose_ );
	pose_.energies().total_energies().show_weighted( std::cerr, score_function_.weights() );
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
		ScoringManager::get_instance()->get_Ramachandran().eval_rama_score_all( pose_, score_function_ );
	}

	if ( check_hbonds ) {
		scoring::hbonds::HBondSet const & hbond_set
			( static_cast< scoring::hbonds::HBondSet const & >
				( pose_.energies().data().get( scoring::EnergiesCacheableDataType::HBOND_SET )));

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


Real
SymAtomTreeMultifunc::operator ()( Multivec const & vars ) const
{
	PROF_START( basic::FUNC );
	symm_min_map_.copy_dofs_to_pose( pose_, vars );
	Real const score( score_function_( pose_ ) );
	PROF_STOP( basic::FUNC );
	return score;
}

void
SymAtomTreeMultifunc::dfunc(
	Multivec const & vars,
	Multivec & dE_dvars
) const {
	PROF_START( basic::DFUNC );

	// in atom_tree_minimize.cc
	atom_tree_dfunc( pose_, symm_min_map_, score_function_, vars, dE_dvars );
	// optional derivative checking
	if ( deriv_check_ ) {
		numerical_derivative_check( symm_min_map_, *this, vars, dE_dvars, deriv_check_verbose_ );
	}
	PROF_STOP( basic::DFUNC );
}

} // namespace symmetry
} // namespace optimization
} // namespace core

