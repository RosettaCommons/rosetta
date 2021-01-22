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
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>
#include <core/scoring/epr_deer/util.hh>

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
static basic::Tracer TR( "core.scoring.epr_deer.metrics.DEERData" );

/// @brief  Print the simulated distance distribution.
/// @param  sim_histr: Simulated distance distribution
/// @param  pose_name: Name of the pose
void
DEERData::print_histogram(
	std::map< Size, Real > const & sim_histr,
	std::string const & pose_name
) const {
	using namespace basic::options;
	if ( option[ OptionKeys::epr_deer::print_data ]() ) {
		for ( auto const & bin_amp : sim_histr ) {
			TR << "\t" << pose_name << "\t";
			for ( auto const & res_sl : residues_ ) {
				TR << res_sl.first << "\t" << res_sl.second << "\t";
			}
			TR << Real( bin_amp.first ) / bins_per_a_ << "\t" <<
				bin_amp.second << std::endl;
		}
	}
}

/// @brief  Function to evaluate score given a distribution
/// @param  sim_histr: Simulated DEER distribution
/// @param  set_score: Boolean to save score or not
/// @return Freshly computed score
Real
DEERData::get_score(
	std::map< Size, Real > const & sim_histr,
	bool const & set_score
) {

	// Get the score (virtual method)
	auto const score = get_score( sim_histr );

	// Set the score if prompted
	if ( set_score ) {
		score_ = score;
	}

	// Return output
	return score;
}

/// @brief  Convolute a distribution with a Gaussian of a specific width
/// @param  distr: DEER Distribution
/// @param  std: Width of Gaussian (mean=0)
/// @return New distribution
std::map< Size, Real >
DEERData::convolute(
	std::map< Size, Real > const & distr,
	Real const & std
) const {

	// Output distribution
	std::map< Size, Real > new_hist;

	// Area under the curve (for normalization)
	Real total = 0.0;

	// Range for iteration
	Size const STDEV_RANGE = 3;
	Real width = STDEV_RANGE * std * bins_per_a_;
	int start = round( distr.begin()->first - width );
	int stop = round( distr.rbegin()->first + width );

	// Outer loop: iterate across distance bins (to add)
	for ( int bin = std::max( 1, start ); bin <= stop + 1; ++bin ) {

		// Calculate amplitude
		Real const fbin = Real( bin ) / bins_per_a_;

		// Inner loop: Iterate across distance bins (to calculate)
		for ( auto const & native_dist : distr ) {

			// Calculate the amplitude
			Real const val_to_add = native_dist.second * gauss(
				fbin, Real( native_dist.first ) / bins_per_a_, std );

			// Add it to the bin
			add_to_map( new_hist, Size( bin ), val_to_add );
			total += val_to_add;
		}
	}

	// Normalize
	for ( auto & new_dist : new_hist ) {
		new_dist.second /= total;
	}

	// Return the distribution
	return new_hist;
}

/// @brief  Returns the residues involved in this data set.
/// @return Vector of residues (ID and spin label type)
/// @detail Residues are saved with two parameters: the residue ID, and
///  the label type. Label type is set to "DEFAULT" by default. Other
///  options include DEFAULT_FAST and CUSTOM
utility::vector1< PairSizeString > const &
DEERData::residues() const {
	return residues_;
}

/// @brief  Returns bins per angstrom for distribution (default: 2)
/// @return Bins per angstrom
Size const &
DEERData::bins_per_a() const {
	return bins_per_a_;
}

/// @brief  Returns the last computed score
/// @return Score (0.0 if never set)
Real
DEERData::score() const {
	return score_;
}

/// @brief  Returns the standard deviation of the distributions to generate
/// @param  Standard deviation (in angstroms)
/// @detail Function has a failsafe to avoid returning a nonzero value
Real
DEERData::stdev() const {
	return std::max( 0.1, stdev_ );
}

/// @brief Sets residue for data set
/// @param residues: Vector of residues (index and spin label type)
void
DEERData::residues(
	utility::vector1< PairSizeString > const & residues
) {
	residues_ = residues;
}

/// @brief Set the number of bins per angstrom for the data set.
/// @param bins_per_a: Bins per angstrom
void
DEERData::bins_per_a(
	Size const & bins_per_a
) {
	bins_per_a_ = bins_per_a;
}

/// @brief Set the score of the data set
/// @param val: Score to save
void
DEERData::score(
	Real const & val
) {
	score_ = val;
}

/// @brief Set the standard deviation of the distributions to generate
/// @param  val: Set the standard deviation to this value
void
DEERData::stdev(
	Real const & val
) {
	stdev_ = val;
}

/// @brief  Computes average distance of distribution for local functions
/// @param  histogram: Distance distribution
/// @return Vector with two items: Average and Standard Deviation
utility::vector1< Real >
DEERData::avg_stdev(
	std::map< Size, Real > const & histogram
) {

	// Initial values for average, variance, and total
	Real bb_avg = 0.0;
	Real bb_var = 0.0;
	Real total = 0.0;

	// Iterate a first time for average
	for ( auto const & p_exp : histogram ) {

		// Aliases
		auto const & dist = p_exp.first / bins_per_a_;
		auto const & amp = p_exp.second;

		// The math
		bb_avg += dist * amp;
		total += amp;
	}

	// Normalize
	bb_avg /= total;

	// Iterate a second time to get variance
	for ( auto const & p_exp : histogram ) {

		// Aliases
		auto const & dist = p_exp.first / bins_per_a_;
		auto const & amp = p_exp.second;

		// The math
		bb_var += pow( ( dist / total ) - bb_avg, 2 ) * amp;
	}

	// Return output
	return utility::vector1< Real >{ bb_avg, sqrt( bb_var ) };
}

/// @brief  Returns the map of distance values used for custom distributions
/// @return Map of values (indeces to distances in Angstroms)
/// @detail Only works when bins_per_a set to zero!
std::map< Size, Real >
DEERData::dist_map() const {
	return dist_map_;
}

/// @brief  Append distance ID to custom distance map
/// @param  dist_id: Unique distance ID used in distance map
/// @param  dist: Distance value in angstroms
/// @detail Only works when bins_per_a set to zero!
void
DEERData::append_dist_id(
	Size dist_id,
	Real dist
) {
	dist_map_[ dist_id ] = dist;
}

}
}
}
}
