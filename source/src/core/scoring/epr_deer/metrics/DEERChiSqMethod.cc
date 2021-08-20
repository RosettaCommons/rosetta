// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERChiSqMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERChiSqMethod.hh>
#include <core/types.hh>



// C++ headers
#include <map>
#include <limits>
#include <cmath> // MANUAL IWYU

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
Real
DEERChiSqMethod::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// The chi-squared function is asymmetric, so the order matters. The
	//  distr_ parameter informs the order
	// ohist: Observed histogram/distribution
	// ehist: Expected histogram/distribution
	auto const & ohist = ( rev_ ) ? sim_histr : distr_;
	auto const & ehist = ( rev_ ) ? distr_ : sim_histr;

	// Set output and iteration range
	Real output = 0.0;

	// Iterate across expected values
	for ( auto const & bin_amp : ehist ) {

		// Define values
		Size const & bin = bin_amp.first;
		auto e_val = std::max( 1e-12, bin_amp.second );
		auto o_val = ( ohist.find( bin ) != ohist.end() )
			? ohist.at( bin ) : std::numeric_limits< Real >::min();

		// Calcualte and add to output
		output += pow( o_val - e_val, 2 ) / e_val;
	}

	// Return
	return sqrt( output / 2.0 );
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
