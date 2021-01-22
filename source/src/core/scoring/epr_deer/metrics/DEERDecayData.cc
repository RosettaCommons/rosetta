// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERData.cc
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Basic headers
#include <basic/datacache/CacheableData.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// C++ headers
#include <string>
#include <set>
#include <map>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.metrics.DEERDecayData" );

/// @brief  Virtual function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @return Freshly computed score
/// @detail A DEER trace is calculated from the simulated histogram.
///  For details, please read del Alamo Biophysical Journal 2020.
///  The score is the average sum-of-squares residuals between the
///  two. Note that there is a slight noise-dependence to this;
///  noisier data will always return higher scores.
Real
DEERDecayData::get_score(
	std::map< Size, Real > const & sim_histr
) {

	// Check if the standard deviation is a free parameter
	if ( fit_stdev_ ) {

		utility::vector1< Real > scores;

		// First get average distance
		std::map< Size, Real > const d
			= { { avg_stdev( sim_histr )[ 1 ], 1 } };

		// Now go through and simulate the DEER traces for each
		for ( Real std = 0.5; std <= 8.0; std += 0.1 ) {
			auto const trace = factory_.trace_from_distr( convolute( d, std ) );
			scores.push_back( sum_of_squares( trace.deer_trace() ) );
		}

		// Return lowest score
		return *std::min_element( scores.begin(), scores.end() );

		// If only the input distribution is used:
	} else {

		// Do the same thing
		auto const trace( factory_.trace_from_distr( sim_histr ) );
		return sum_of_squares( trace.deer_trace() );
	}
}

/// @brief Calculate sum of squares of simulated DEER trace
/// @param  trace1: Simulated DEER trace
/// @param  trace2: Reference DEER trace; experimental by default
/// @param  normalize: Normalize by number of time points; recommended
/// @return Sum of squared residuals
Real
DEERDecayData::sum_of_squares(
	utility::vector1< Real > const & sim_trace,
	bool const & normalize // = true
) const {

	// Define output;
	Real residuals = 0.0;

	// Define trace for comparison (trace2 size is zero if not passed)
	auto const & exp_trace = factory_.trace();

	// Make sure everything is in order
	debug_assert( sim_trace.size() == exp_trace.size() );

	// Iterate over traces and get residuals
	for ( Size i = 1; i <= sim_trace.size(); ++i ) {
		residuals += pow( ( sim_trace[ i ] - exp_trace[ i ] ) / noise_, 2 );
	}

	// Return mean sum of squares (or sum of squares if normalize == false)
	return ( normalize ) ? residuals / sim_trace.size() : residuals;
}

/// @brief Initialize DEER factory object
/// @param trace: Experimental DEER trace
/// @param tiem_pts: Time points corresponding to DEER trace
void
DEERDecayData::init_factory(
	utility::vector1< Real > const & trace,
	utility::vector1< Real > const & time_pts
) {
	factory_ = Simulated4PDEERTraceFactory(
		trace, time_pts, bckg_, bins_per_a_, max_dist_ );
}

/// @brief Returns DEER trace factory object
/// @return Factory object
Simulated4PDEERTraceFactory
DEERDecayData::factory() const{
	return factory_;
}

/// @brief  Returns the noise from the imaginary component
/// @return Noise level
Real const &
DEERDecayData::noise() const {
	return noise_;
}

/// @brief  Returns whether standard deviation is a fitting parameter
/// @return Boolean
bool const &
DEERDecayData::fit_stdev() const {
	return fit_stdev_;
}

/// @brief Returns background type
/// @return Intermolecular coupling background type
std::string
DEERDecayData::bckg() const {
	return bckg_;
}

/// @brief Returns maximum distance for kernel calculation
/// @return Maximum distance for kernel calculation
Size
DEERDecayData::max_dist() const {
	return max_dist_;
}

/// @brief  Sets DEER trace factory object
/// @param  factory: Factory object
void
DEERDecayData::factory(
	Simulated4PDEERTraceFactory const & val
) {
	factory_ = val;
}

/// @brief Sets the noise from the imaginary component
/// @param val: Noise level
void
DEERDecayData::noise(
	Real const & val
) {
	noise_ = val;
}

/// @brief Sets if standard deviation is varied as a parameter
/// @param val: Boolean
void
DEERDecayData::fit_stdev(
	bool const & val
) {
	fit_stdev_ = val;
}

/// @brief Sets background type
/// @param  val: Type of intermolecular coupling background
void
DEERDecayData::bckg(
	std::string const & val
) {
	bckg_ = val;
}

/// @brief Sets maximum distance for kernel calculation
/// @param val: Maximum distance to set
void
DEERDecayData::max_dist(
	Size const & val
) {
	max_dist_ = val;
}

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core
