// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERJSMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.hh>
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
DEERJSMethod::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// Initialize bins and score
	Real jensen_shannon_d = 0.0;
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

		// Set the middle distance and add to the score
		Real const mid = ( sim + exp ) / 2.0;
		jensen_shannon_d += sim * ln( sim / mid ) + exp * ln( exp / mid );
	}

	// Returnt he score
	return sqrt( jensen_shannon_d / 2.0 );
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
