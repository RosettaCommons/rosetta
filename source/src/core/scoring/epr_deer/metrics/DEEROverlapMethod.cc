// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEEROverlapMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEEROverlapMethod.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/util.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

#include <boost/math/special_functions/erf.hpp>

// C++ headers
#include <string>
#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
Real
DEEROverlapMethod::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// Initialize some variables first
	Real score = 0.0;
	Real supremum = 0.0;
	Real exp_integral_lo = 0.0;
	Real exp_integral = 0.0;
	Real exp_integral_hi = 0.0;
	Real sim_integral = 0.0;

	// Set minimum and maximum bins
	Size min_bin = 0;
	Size max_bin = 0;
	std::tie( min_bin, max_bin ) = min_max( sim_histr, distr_ );

	// Iterate across bins
	for ( Size bin = min_bin; bin <= max_bin; ++bin ) {

		// Identify amplitudes at these bins
		bool const sim = sim_histr.find( bin ) != sim_histr.end();
		bool const exp = distr_.find( bin ) != distr_.end();

		// If integrals must be computed (for example, for the Wasserstein
		// and Kolmogorov-Smirnov metrics)
		if ( integral_ ) {
			sim_integral  += ( sim ) ? sim_histr.at( bin )  : 0.0;
			exp_integral_lo += ( exp ) ? lower_bound_.at( bin ) : 0.0;
			exp_integral  += ( exp ) ? distr_.at( bin )   : 0.0;
			exp_integral_hi += ( exp ) ? upper_bound_.at( bin ) : 0.0;

			// If confidence bands are not being used
			if ( !bounds_ ) {
				auto const val = std::abs( sim_integral - exp_integral );

				// Score is simply the area between these values
				// Supremum is the maximum such discrepancy
				score += val;
				if ( val > supremum ) {
					supremum = val;
				}
			} else {

				// Score is the extent to which this exceeds confidence bands
				auto const val_upper = std::max( 0.0, sim_integral - exp_integral_hi );
				auto const val_lower = std::max( 0.0, exp_integral_lo - sim_integral );
				score += val_lower + val_upper;
				if ( val_upper + val_lower > supremum ) {
					supremum = val_upper + val_lower;
				}
			}

			// If integrals are not being used
		} else {

			// Fetch values and compare
			auto const sim_val = ( sim ) ? sim_histr.at( bin ) : 0.0;
			auto const exp_val = ( exp ) ? ( ( bounds_ )
				? upper_bound_.at( bin ) : distr_.at( bin ) ) : 0.0;

			// Update supremum / greatest discrepancy
			if ( std::abs( sim_val - exp_val ) > supremum ) {
				supremum = std::abs( sim_val - exp_val );
			}

			// Deduct overlap from score
			// Must be deducted since we want to reward overlap
			score -= ( sim && exp )? std::min( sim_val, exp_val ) : 0.0;
		}
	}

	// For Kolmogorov-Smirnov or Maximum Discrepancy
	if ( singleval_ ) {
		return supremum;

		// For Overlap or Wasserstein
	} else {
		return score;
	}
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
