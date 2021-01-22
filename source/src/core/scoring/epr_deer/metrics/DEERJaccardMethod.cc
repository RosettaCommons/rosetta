// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERJaccardMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERJaccardMethod.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/util.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
Real
DEERJaccardMethod::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// Initialize bins
	Real numerator = 0.0;
	Real denominator = 0.0;
	Size min_bin = 0;
	Size max_bin = 0;
	std::tie( min_bin, max_bin ) = min_max( sim_histr, distr_ );

	// Iterate through bins
	for ( Size bin = min_bin; bin <= max_bin; ++bin ) {

		Real const & sim = ( sim_histr.find( bin ) != sim_histr.end() )
			? sim_histr.at( bin ) : std::numeric_limits< Real >::min();

		// If confidence bands are being used
		if ( bounds_ ) {
			if ( upper_bound_.find( bin ) != upper_bound_.end() ) {

				// Compare to lower and upper bounds
				numerator += std::min( sim, upper_bound_.at( bin ) );
				denominator += std::max( sim, lower_bound_.at( bin ) );
			} else {
				denominator += sim;
			}

			// Otherwise, just compare to best fit
		} else {
			if ( distr_.find( bin ) != distr_.end() ) {
				numerator += std::min( sim, distr_.at( bin ) );
				denominator += std::max( sim, distr_.at( bin ) );
			} else {
				denominator += sim;
			}
		}
	}

	// Return a negative score, since we want to maximize this value
	return -1 * ( numerator / denominator );
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
