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
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>

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

static basic::Tracer TR( "core.scoring.epr_deer.DEERData" );

/// @brief Print the simulated distance distribution.
void
DEERData::print_histogram(
	std::map< Size, Real > const & simulated_histogram
) const {
	if ( basic::options::option[ basic::options::OptionKeys::epr_deer::print_data ]() ) {
		TR.Trace << "For dataset containing the following residues:" << std::endl;
		for ( auto const & res_sl : residues_ ) {
			TR.Trace << "\t" << res_sl.first << "\t" << res_sl.second << std::endl;
		}
		TR.Trace << "Distance distribution:" << std::endl;
		for ( auto const & bin_amp : simulated_histogram ) {
			TR.Trace << "\t" << Real( bin_amp.first ) / bins_per_angstrom_ << "\t" << bin_amp.second << std::endl;
		}
		TR.Trace << std::endl;
	}
}

/// @brief  Virtual scoring method. Returns value of zero for the base class.
/// @detail Takes as input a simulated histogram and an option to set the private
///     score value to the result (default: false. This saves some code)
Real
DEERData::get_score(
	std::map< Size, Real > const &,
	bool const & set_score
) {
	if ( set_score ) {
		score_ = 0.0;
	}
	return 0.0;
}

/// @brief Returns the residues involved in this data set.
/// @detail Residues are saves with two parameters: the residue ID, and the label type
///     Label type is set to "DEFAULT" by default (duh). Other options include
///     DEFAULT_FAST and CUSTOM
utility::vector1< std::pair< Size, std::string > > const &
DEERData::residues() const {
	return residues_;
}

/// @brief Returns the granularity of the distance distribution. Default is 2 bins per A
Size const &
DEERData::bins_per_angstrom() const {
	return bins_per_angstrom_;
}

/// @brief Returns the last computed score. Obtained using get_score() or manually set
Real
DEERData::score() const {
	return score_;
}

/// @brief Returns the relative weight. Default: 1. Can be lowered when the data is less trustworthy
Real
DEERData::relative_weight() const {
	return rel_weight_;
}
/*
/// @brief Returns to condition ID to which this data belongs
Size
DEERData::condition() const {
return condition_;
}
*/
/// @brief Sets residue for data set; info for each "residue" consists of the index and the spin label type
void
DEERData::residues(
	utility::vector1< std::pair< Size, std::string > > const & residues
) {
	residues_ = residues;
}

/// @brief Set the number of bins per angstrom for the data set.
void
DEERData::bins_per_angstrom(
	Size const & bins_per_angstrom
) {
	bins_per_angstrom_ = bins_per_angstrom;
}

/// @brief Set the score of the data set
void
DEERData::score(
	Real const & val
) {
	score_ = val;
}

/// @brief Set the relative weight of the data set
void
DEERData::relative_weight(
	Real const & val
) {
	rel_weight_ = val;
}
/*
/// @brief Returns to condition ID to which this data belongs
void
DEERData::condition(
Size const & val
) {
condition_ = val;
}
*/
/// @brief Computes average distance of distribution for local functions
Real
DEERData::compute_avg_dist(
	std::map< Size, Real > const & simulated_histogram
) const {
	Real total_dist = 0.0;
	Real total_auc  = 0.0;
	for ( auto const & dist : simulated_histogram ) {
		total_dist += dist.first * dist.second / Real( bins_per_angstrom_ );
		total_auc  += dist.second;
	}
	return total_dist / total_auc;
}

//////////////////////////////////

/// @brief Computes average distance and score based on bounds.
Real
DEERDistanceBounds::get_score(
	std::map< Size, Real > const & simulated_histogram,
	bool const & set_score
) {
	Real average = compute_avg_dist( simulated_histogram );
	Real score_lower = std::max( bounds_.first - average,  0.0 ) / step_;
	Real score_upper = std::max( average - bounds_.second, 0.0 ) / step_;
	if ( set_score ) {
		score_ = pow( score_lower + score_upper, 2 );
	}
	return pow( score_lower + score_upper, 2 );
}

/// @brief Returns the lower and upper distance bounds
std::pair< Real, Real > const &
DEERDistanceBounds::bounds() const {
	return bounds_;
}

/// @brief Returns the step / steepness
Real const &
DEERDistanceBounds::step() const {
	return step_;
}

/// @brief Sets the lower and upper bounds
void
DEERDistanceBounds::bounds(
	Real const & lower,
	Real const & upper
) {
	auto const & lower_new = std::min( lower, upper );
	auto const & upper_new = std::max( lower, upper );
	bounds_ = std::make_pair( lower_new, upper_new );
}

/// @brief Sets the step
void
DEERDistanceBounds::step(
	Real const & step
) {
	step_ = step;
}

//////////////////////////////////

/// @brief  Computes cross-entropy of simulated distribution from experimental.
Real
DEERDistanceDistribution::get_score(
	std::map< Size, Real > const & simulated_histogram,
	bool const & set_score
) {
	Real score = 0.0;
	for ( auto const & dist_vals : simulated_histogram ) {
		Real const & dist = dist_vals.first;
		Real const & val = dist_vals.second;
		Real exp_val = ( best_fit_.find( dist ) == best_fit_.end() )
			? std::abs( std::numeric_limits< Real >::min() ) : best_fit_.at( dist );
		// Here to avoid getting nan or inf as a score
		if ( std::isinf( std::abs( 1.0 / exp_val ) ) ) {
			score += val * log( std::abs( std::numeric_limits< Real >::max() ) );
		} else {
			score += val * log( 1.0 / exp_val );
		}
	}
	if ( set_score ) {
		score_ = score;
	}
	return score;
}

/// @brief Returns the lower bound/confidence band for the distance distribution
std::map< Size, Real > const &
DEERDistanceDistribution::lower_bound() const {
	return lower_bound_;
}

/// @brief Returns the line of best fit. Used to calculate cross-entropy
std::map< Size, Real > const &
DEERDistanceDistribution::best_fit() const {
	return best_fit_;
}

/// @brief Returns the upper bound/confidence band for the distance distribution
std::map< Size, Real > const &
DEERDistanceDistribution::upper_bound() const {
	return upper_bound_;
}

/// @brief Sets the lower bound/confidence band for the distance distribution
void
DEERDistanceDistribution::lower_bound(
	std::map< Size, Real > const & lower_bound
) {
	lower_bound_ = lower_bound;
}

/// @brief Sets the line of best fit. Used to calculate cross-entropy
void
DEERDistanceDistribution::best_fit(
	std::map< Size, Real > const & val
) {
	best_fit_ = val;
}

/// @brief Sets the upper bound/confidence band for the distance distribution
void
DEERDistanceDistribution::upper_bound(
	std::map< Size, Real > const & upper_bound
) {
	upper_bound_ = upper_bound;
}

/// @brief Adds data to map at a specific distance. Overwrites if distance is occupied
void
DEERDistanceDistribution::append(
	Size const & dist_bin,
	Real const & lower,
	Real const & upper
) {
	lower_bound_[ dist_bin ] = std::min( lower, upper );
	upper_bound_[ dist_bin ] = std::max( lower, upper );
}

//////////////////////////////////

/// @brief Computes sum-of-squares of experimental to simulated DEER trace
Real
DEERDecayData::get_score(
	std::map< Size, Real > const & simulated_histogram,
	bool const & set_score
) {
	Real score = 0.0;
	Simulated4PDEERTrace sim_trace;
	// If the standard deviation of the distribution is fitted as a parameter
	if ( fit_stdev_ ) {
		Real best_score = 0.0;
		Real avg = compute_avg_dist( simulated_histogram );
		for ( Real stdev = 1.0; stdev <= 8.0; stdev += 0.5 ) {
			auto sim_histogram_mod = mod_distr( avg, stdev );
			sim_trace.simulate_decay( trace_, fit_info_, sim_histogram_mod, bins_per_angstrom_ );
			auto sum_of_squares = sim_trace.sum_of_squares( trace_, noise_ );
			if ( best_score == 0.0 || sum_of_squares < best_score ) {
				best_score = sum_of_squares;
			}
		}
		score = best_score;
		// If the entire distribution is being used
	} else {
		sim_trace.simulate_decay( trace_, fit_info_, simulated_histogram, bins_per_angstrom_ );
		score = sim_trace.sum_of_squares( trace_, noise_ );
	}
	if ( set_score ) {
		score_ = score;
	}
	// Print the DEER trace if requested (default: false)
	if ( basic::options::option[ basic::options::OptionKeys::epr_deer::print_data ]() ) {
		TR.Trace << "DEER trace:" << std::endl;
		for ( auto const & time_trace : trace_ ) {
			TR.Trace << "\t" << time_trace.first << "\t" << sim_trace[ time_trace.first ] << "\t" << time_trace.second << std::endl;
		}
		TR.Trace << std::endl;
	}
	return score;
}

/// @brief Returns a gaussian distribution with certain avg and stdev.
std::map< Size, Real >
DEERDecayData::mod_distr(
	Real const & avg,
	Real const & stdev
) const {
	Real N_STDEVS = 4.0; // Number of standard deviations - 4 is overkill.
	std::map< Size, Real > output;
	Size min_dist = round( std::max( 1.0, avg - ( N_STDEVS * stdev ) ) ) * bins_per_angstrom_;
	Size max_dist = round( avg + ( N_STDEVS * stdev ) ) * bins_per_angstrom_;
	Real baseline = 0.0;
	for ( Size dist = min_dist; min_dist <= max_dist; ++min_dist ) {
		Real mod_dist = Real( dist ) / bins_per_angstrom_;
		// Gaussian distribution
		output[ dist ] = std::exp( -0.5 * ( pow( ( mod_dist - avg ) / stdev, 2 ) ) ) /
			( stdev * pow( 2.0 * numeric::constants::d::pi, 0.5 ) );
		baseline += output[ dist ];
	}
	for ( auto & dist : output ) {
		dist.second /= baseline;
	}
	return output;
}

/// @brief Returns the experimental DEER trace. Key is time point in microseconds.
std::map< Real, Real > const &
DEERDecayData::trace() const {
	return trace_;
}

/// @brief Returns the last value obtained when fitting background slope (stored in fitting info)
Real const &
DEERDecayData::k_fit() const {
	return fit_info_.last_slope_;
}

/// @brief Returns the last value obtained when fitting modulation depth (stored in fitting info)
Real const &
DEERDecayData::modulation_depth_fit() const {
	return fit_info_.last_mod_depth_;
}

/// @brief Returns the time points squared, used to fit slope and avoids unnecessary calculation.
Real const &
DEERDecayData::time_points_sqd() const {
	return fit_info_.time_pts_sqd_;
}

/// @brief Returns the spin value for a given time point and distance. Computed at the front end.
Real const &
DEERDecayData::spin_val(
	Size const & dist_bin,
	Real const & time_pt
) const {
	return fit_info_.spin_map_.at( dist_bin ).at( time_pt );
}

/// @brief Returns the upper and lower bounds of the modulation depth
std::pair< Real, Real > const &
DEERDecayData::mod_depth_bounds() const {
	return fit_info_.mod_depth_bounds_;
}

std::string const &
DEERDecayData::bckg() const {
	return fit_info_.bckg_;
}

/// @brief Returns the noise from the imaginary component, provided manually as an option
Real const &
DEERDecayData::noise() const {
	return noise_;
}

/// @brief Returns whether or not the standard deviation is used as a fitting parameter
bool const &
DEERDecayData::fit_stdev() const {
	return fit_stdev_;
}

/// @brief Sets the experimental DEER trace
void
DEERDecayData::trace(
	std::map< Real, Real > const & trace
) {
	trace_ = trace;
}

/// @brief Appends time point and signal values for the experimental DEER trace
void
DEERDecayData::append_trace_data(
	Real const & time,
	Real const & signal
) {
	trace_[ time ] = signal;
}

/// @brief Sets the value obtained when fitting background slope (stored in fitting info)
void
DEERDecayData::k_fit(
	Real const & k
) {
	fit_info_.last_slope_ = k;
}

/// @brief Sets the value obtained when fitting modulation depth (stored in fitting info)
void
DEERDecayData::modulation_depth_fit(
	Real const & modulation_depth
) {
	fit_info_.last_mod_depth_ = modulation_depth;
}

/// @brief Sets the value for time points squared, used to fit slope and avoids unnecessary calculation.
void
DEERDecayData::time_pts_sqd(
	Real const & time_pts_sqd
) {
	fit_info_.time_pts_sqd_ = time_pts_sqd;
}

/// @brief Sets the "spin map" that contains trace values at all time points and distances
void
DEERDecayData::spin_map(
	std::map< Size, std::map< Real, Real > > const & spin_map
) {
	fit_info_.spin_map_ = spin_map;
}

/// @brief Sets lower and upper bounds for modulation depth
void
DEERDecayData::mod_depth_bounds(
	Real const & lower,
	Real const & upper
) {
	auto const & lower_new = std::min( lower, upper );
	auto const & upper_new = std::max( lower, upper );
	fit_info_.mod_depth_bounds_ = std::make_pair( lower_new, upper_new );
}

/// @brief Sets what type of background is used for fitting: options are 3D, NON-3D, NONE
void
DEERDecayData::bckg(
	std::string const & bckg
) {
	fit_info_.bckg_ = bckg;
}

/// @brief Sets the noise from the imaginary component, provided manually as an option
void
DEERDecayData::noise(
	Real const & noise
) {
	noise_ = noise;
}

/// @brief Sets whether the standard deviation is a parameter when fitting DEER traces
void
DEERDecayData::fit_stdev(
	bool const & fit_stdev
) {
	fit_stdev_ = fit_stdev;
}

/// @brief Appends data to experimental decay trace and values to relevant variables
void
DEERDecayData::append_trace_data_and_calculate(
	Real const & time,
	Real const & signal,
	Size const & max_bin
) {
	// First append the data to the DEER trace
	trace_[ time ] = signal;

	// Then add the time points squared to the appropriate private data member
	fit_info_.time_pts_sqd_ += pow( time, 2 );

	// Finally add the data to the spin map
	Size min_bin = round( 15.0 * bins_per_angstrom_ );
	for ( Size bin = min_bin; bin <= max_bin; ++bin ) {
		Real r = Real( bin ) / bins_per_angstrom_;
		if ( fit_info_.spin_map_.find( bin ) == fit_info_.spin_map_.end() ) {
			fit_info_.spin_map_[ bin ] = std::map< Real, Real >();
		}
		fit_info_.spin_map_[ bin ][ time ] = calculate_decay_2_spin( time, r );
	}
	if ( fit_info_.spin_map_[ min_bin ].find( 0.0 ) == fit_info_.spin_map_[ min_bin ].end() ) {
		append_trace_data_and_calculate( 0.0, 1.0, max_bin );
	}
}

/// @brief Calculates the DEER trace given a time point (microseconds) and distance (A)
Real
DEERDecayData::calculate_decay_2_spin(
	Real const & t,
	Real const & r
) const {
	Size MAX_BINS = 201;
	Real output( 0.0 );
	// Compute the intramolecular DEER signal.
	// The variable ii is the angle between the interspin vector and the external magnetic field
	for ( Size ii = 0 ; ii < MAX_BINS; ++ii ) {
		Real theta = Real( ii * numeric::constants::d::pi_over_2 ) / Real( 201 - 1 );
		Real term1 = 326.98 * ( 1.0 - 3.0 * pow( std::cos( theta ), 2 ) );
		Real term2 = std::sin( theta ) * std::cos( ( term1 * t ) / pow( r / 10.0, 3 ) );
		output += term2;
	}
	return ( output / MAX_BINS );
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
