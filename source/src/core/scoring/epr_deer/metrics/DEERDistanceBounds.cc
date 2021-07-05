// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/DEERData.cc
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.hh>


// Project headers
#include <core/types.hh>

// Basic headers

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// Basic headers

// C++ headers
#include <map>

#include <cmath> // MANUAL IWYU

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
Real
DEERDistanceBounds::get_score(
	std::map< Size, Real > const & sim_histr
) {
	Real const average = avg_stdev( sim_histr )[ 1 ];
	Real const lower = std::max( lo_ - average,  0.0 ) / step_;
	Real const upper = std::max( average - hi_, 0.0 ) / step_;
	return pow( lower + upper, 2 );
}

/// @brief  Returns the lower and upper distance bounds
/// @return Pair of lower and upper bounds
std::pair< Real, Real >
DEERDistanceBounds::bounds() const {
	return std::make_pair( lo_, hi_ );
}

/// @brief  Returns the step / steepness of the scoring function
/// @return Steepness value
Real const &
DEERDistanceBounds::step() const {
	return step_;
}

/// @brief Sets the lower and upper bounds
/// @param lo: Lower bound
/// @param hi: Upper bound
void
DEERDistanceBounds::bounds(
	Real const & lo,
	Real const & hi
) {
	lo_ = std::min( lo, hi );
	hi_ = std::max( lo, hi );
}

/// @brief Sets the step / steepness
/// @param step: Step
void
DEERDistanceBounds::step(
	Real const & step
) {
	step_ = step;
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
