// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERMiscMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.hh>
#include <core/scoring/epr_deer/util.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>



// C++ headers
#include <string>
#include <map>
#include <limits>
#include <tuple>
#include <cmath> // MANUAL IWYU

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
Real
DEERMiscMethod::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// Initialize bins
	Real score = 0.0;
	Size min_bin = 0;
	Size max_bin = 0;
	std::tie( min_bin, max_bin ) = min_max( sim_histr, distr_ );

	// Iterate through bins
	for ( Size bin = min_bin; bin <= max_bin; ++bin ) {

		// If these bins don't exist in either distribution, replace
		// them with an arbitrarily low value
		Real const & sim = ( sim_histr.find( bin ) != sim_histr.end() )
			? sim_histr.at( bin ) : std::numeric_limits< Real >::min();

		Real const & exp = ( distr_.find( bin ) != distr_.end() )
			? distr_.at( bin ) : std::numeric_limits< Real >::min();

		// Add to score, depending on scoring method
		score += ( mode_ == "CONDITIONAL" ) ? exp * sim : sqrt( exp * sim );
	}

	// Now return some function of the score, depending on scoring method
	score = std::max( score, std::numeric_limits< Real >::min() );
	return ( mode_ == "HELLINGER" ) ? 1.0 - score : -1 * log( score );

}

/// @brief Set mode being used to score
/// @param  val: Mode with which to score
void
DEERMiscMethod::mode(
	std::string const & val
) {
	if ( std::none_of( modes_.begin(), modes_.end(),
			[&]( std::string const & s ){ return s == val; } )
			) {
		throw CREATE_EXCEPTION( utility::excn::KeyError,
			"Mode type " + val + " not found!" );
	} else {
		mode_ = val;
	}
}

/// @brief Return the mode being used to score
/// @return Mode being used to score
std::string
DEERMiscMethod::mode() const {
	return mode_;
}


} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
