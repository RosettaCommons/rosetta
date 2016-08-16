// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMLJScore.cc
/// @brief  Molecular mechanics lj score class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMLJScore.hh>
#include <core/scoring/mm/MMLJLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/scoring/ScoringManager.hh>


// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>
#include <math.h>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace mm {

MMLJScore::MMLJScore() :
	mm_lj_library_( scoring::ScoringManager::get_instance()->get_MMLJLibrary() )
{ }

MMLJScore::MMLJScore( MMLJLibrary const & mmljl ) :
	mm_lj_library_( mmljl )
{ }

MMLJScore::~MMLJScore() {}

/// @details blah
Energy
MMLJScore::score( Size atom1, Size atom2, Size path_distance, Real distance_sq ) const
{
	// lookup params
	mm_lj_param_set atom1_params, atom2_params;
	if ( path_distance == 3 ) {
		atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
		atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
	} else {
		atom1_params = mm_lj_library_.lookup( atom1 );
		atom2_params = mm_lj_library_.lookup( atom2 );
	}

	// calc score
	Real epsilon( sqrt( atom1_params.key2() * atom2_params.key2() ) );
	Real rminsq_over_distsq( ( atom1_params.key1() + atom2_params.key1() ) * ( atom1_params.key1() + atom2_params.key1() ) / distance_sq );
	Real rminsq_over_distsq_third = rminsq_over_distsq * rminsq_over_distsq * rminsq_over_distsq;

	Real score = epsilon * rminsq_over_distsq_third * ( rminsq_over_distsq_third - 2 );

	return score;
}

/// @details blah
Energy
MMLJScore::deriv_score( Size atom1, Size atom2, Size path_distance, Real distance_sq ) const
{
	// lookup params
	mm_lj_param_set atom1_params, atom2_params;
	if ( path_distance == 3 ) {
		atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
		atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
	} else {
		atom1_params = mm_lj_library_.lookup( atom1 );
		atom2_params = mm_lj_library_.lookup( atom2 );
	}

	// calc deriv
	Real epsilon( sqrt( atom1_params.key2() * atom2_params.key2() ) );
	Real rminsq( ( atom1_params.key1() + atom2_params.key1() ) * ( atom1_params.key1() + atom2_params.key1() ) );
	Real rmin_sixth = rminsq * rminsq * rminsq;
	Real rmin_twelfth = rmin_sixth * rmin_sixth;
	Real inv_dist_sq = 1/distance_sq;
	Real inv_dist_sixth = inv_dist_sq * inv_dist_sq * inv_dist_sq;
	Real deriv = 6 * epsilon * inv_dist_sixth * ( rmin_sixth * inv_dist_sq - rmin_twelfth * inv_dist_sixth * inv_dist_sq );

	return deriv;
}

/// @details Notably, this returns the minimum distance, NOT the minimum dist sq
Real
MMLJScore::min_dist( Size atom1, Size atom2, Size path_distance ) const
{
	// lookup params
	mm_lj_param_set atom1_params, atom2_params;
	if ( path_distance == 3 ) {
		atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
		atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
	} else {
		atom1_params = mm_lj_library_.lookup( atom1 );
		atom2_params = mm_lj_library_.lookup( atom2 );
	}

	// calc min
	Real rmin ( atom1_params.key1() + atom2_params.key1() );

	return rmin;
}

} // namespace mm
} // namespace scoring
} // namespace core
