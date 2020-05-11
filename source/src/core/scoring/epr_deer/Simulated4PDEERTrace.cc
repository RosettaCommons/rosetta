// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/decay/Simulated4PDEERTrace.cc
/// @brief  Class that simulates the observable signal resulting from electron dipolar coupling
/// @details The dipolar coupling signal between two electrons can be simulated from a distance
///       distribution. This class does that simulation by pulling values from the datacache
///       container objects (DEERData) and effectively storing a vector with values of interest.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>

#include <core/scoring/epr_deer/DEERData.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <map>

namespace core {
namespace scoring {
namespace epr_deer {

static basic::Tracer TR( "core.scoring.epr_deer.Simulated4PDEERTrace" );

/// @bried  Default constructor
Simulated4PDEERTrace::Simulated4PDEERTrace() {}

/// @brief Custom constructor
Simulated4PDEERTrace::Simulated4PDEERTrace(
	std::map< Real, Real > const & exp_trace,
	FittingInfo & fit_info,
	std::map< Size, Real > const & sim_histogram,
	Size const & bins_per_a
) {
	simulate_decay( exp_trace, fit_info, sim_histogram, bins_per_a );
}

/// @brief Bracket operator for accessing the DEER trace
Real
Simulated4PDEERTrace::operator[](
	Real const & t
) const {
	if ( trace_.find( t ) == trace_.end() ) {
		return 0.0;
	} else {
		return trace_.at( t );
	}
}

/// @brief Main function for simulating 4-pulse DEER trace
/// @detail This function first simulates the intramolecular dipolar coupling
///     form factor using the evolution kernel described in
///     "DEER Distance Measurements on Proteins" by Jeschke (2012).
///     It then finds the best background fit for the experimental data,
///     if necessary (3D by default, but can be set to "NONE" or "NON-3D")
///     Background-fitted DEER trace is saved and returned.
std::map< Real, Real >
Simulated4PDEERTrace::simulate_decay(
	std::map< Real, Real > const & exp_trace,
	FittingInfo & fit_info,
	std::map< Size, Real > const & sim_histogram,
	Size const & bins_per_a
) {
	auto trace_no_bckg = calc_4pdeer( exp_trace, fit_info, sim_histogram, bins_per_a );
	trace_ = optimize_bckg( exp_trace, trace_no_bckg, fit_info );
	return trace_;
}

/// @brief Deterine the intramolecular DEER trace
std::map< Real, Real >
Simulated4PDEERTrace::calc_4pdeer(
	std::map< Real, Real > const & exp_trace,
	FittingInfo const & fit_info,
	std::map< Size, Real > const & sim_histogram,
	Size const & bins_per_a
) {
	std::map< Real, Real > output;
	// No DEER signal is possible for distances closer than 15 angstroms due to short-distance coupling
	Size min_dist = std::max( 15 * bins_per_a, sim_histogram.begin()->first );
	Size max_dist = std::min( fit_info.spin_map_.rbegin()->first, sim_histogram.rbegin()->first );
	// Initialize the array
	for ( auto const & time_signal : exp_trace ) {
		output[ time_signal.first ] = 0.0;
	}
	// Outer loop: iterate through each distance
	// This order reduces the number of times the if statement needs to be checked
	for ( Size r = min_dist; r <= max_dist; ++r ) {
		if ( sim_histogram.find( r ) == sim_histogram.end() ) {
			continue;
		}
		// inner loop: iterate through each time point
		for ( auto const & time_signal : exp_trace ) {
			output[ time_signal.first ] += fit_info.spin_map_.at( r ).at( time_signal.first );
		}
	}
	// Now normalize for time point zero
	Real baseline = output[ 0.0 ];
	if ( baseline == 0.0 ) {
		output = fit_info.spin_map_.rbegin()->second;
	}
	for ( auto & time_signal : output ) {
		time_signal.second /= baseline;
	}
	return output;
}

/// @brief Optimize the background function of the DEER trace: general function
std::map< Real, Real >
Simulated4PDEERTrace::optimize_bckg(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & trace_no_bckg,
	FittingInfo & fit_info
) {
	if ( fit_info.bckg_ == "NONE" ) {
		return trace_no_bckg;
	} else if ( fit_info.bckg_ == "3D" ) {
		return optimize_bckg( exp_trace, trace_no_bckg, fit_info, false );
	} else if ( fit_info.bckg_ == "NON_3D" ) {
		return optimize_bckg( exp_trace, trace_no_bckg, fit_info, true );
	} else {
		utility_exit_with_message( "DEER BACKGROUND FXN NOT DECLARED! QUITTING." );
	}
}

/// @brief Fit the background decay using linear regression
Real
Simulated4PDEERTrace::fit_k(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & trace,
	FittingInfo & fit_info,
	Real const & d,
	Real const & l
) {
	// This function assumes the y-intercept equals zero. Theoretically this must be the case
	Real xy = 0.0;
	Real xx = ( d == 3.0 ) ? fit_info.time_pts_sqd_ : 0.0;
	for ( auto const & time_exp : exp_trace ) {
		Real sim = ( 1.0 - l * ( 1.0 - trace.at( time_exp.first ) ) );
		xy += pow( std::abs( time_exp.first ), ( d / 3.0 ) )
			* log( time_exp.second / sim );
		if ( d != 3.0 ) {
			xx += pow( pow( std::abs( time_exp.first ), ( d / 3.0 ) ), 2.0 );
		}
	}
	return xy / xx;
}
/// @brief Optimize the modulation depth using simple gradient descent
void
Simulated4PDEERTrace::optimize_mod_depth(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & trace,
	FittingInfo & fit_info,
	Real & d,
	Real & k,
	Real & l,
	bool fit_dim // = false
) {
	Real END_STEP = 0.00025;
	Size n_steps = 0;
	Real step = 0.0;
	do {
		++n_steps;
		step = 0.0;
		k = fit_k( exp_trace, trace, fit_info, d, l );
		for ( auto const & time_exp : exp_trace ) {
			// This is the derivative of the background function with respect to l
			Real bckg = k * pow( std::abs( time_exp.first ), ( d / 3.0 ) );
			Real term1 = -1.0 * std::exp( bckg );
			Real term2 = std::exp( bckg ) * ( 1.0 - l * ( 1.0 - trace.at( time_exp.first ) ) ) - time_exp.second;
			Real term3 = -1.0 * trace.at( time_exp.first ) + 1.0;
			step += term1 * term2 * term3;
		}
		l -= step / trace.size();
		if ( fit_dim && d > 0.01 ) {
			// The analytical derivative of this function with respect to d has the term ln( kt )
			// Since kt is always negative we will do this incrementally.
			k = fit_k( exp_trace, trace, fit_info, d, l );
			auto temp_trace_minus = fit_trace( exp_trace, trace, d - 0.01, k, l );
			Real sum_of_squares_minus = sum_of_squares( exp_trace, temp_trace_minus );
			auto temp_trace_plus = fit_trace( exp_trace, trace, d + 0.01, k, l );
			Real sum_of_squares_plus = sum_of_squares( exp_trace, temp_trace_plus );
			if ( sum_of_squares_minus < sum_of_squares_plus && !std::isnan( std::abs( sum_of_squares_minus ) ) ) {
				d -= 0.01;
			} else if ( sum_of_squares_minus > sum_of_squares_plus && !std::isnan( std::abs( sum_of_squares_plus ) ) ) {
				d += 0.01;
			}
		}
	} while ( l >= fit_info.mod_depth_bounds_.first
			&& l <= fit_info.mod_depth_bounds_.second
			&& step / trace.size() > END_STEP
			&& n_steps < 100
			&& d <= 3.5
			&& d >= 2.0
			);
}

/// @brief Optimize the background function of the DEER trace: detailed function
// @detail  This is a nonlinear least-squares "fit-within-a-fit" approach for
//     determining the background parameters. To save time we save the results
//     after solving them. If these parameters have not previously been optimized,
//     an initial search is performed (the first if loop)
std::map< Real, Real >
Simulated4PDEERTrace::optimize_bckg(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & trace,
	FittingInfo & fit_info,
	bool const & optimize_dim
) {

	Real k = fit_info.last_slope_;
	Real l = fit_info.last_mod_depth_;
	Real d = fit_info.last_dim_;
	// Perform initial brute-force search if the parameters were not initialized.
	if ( k == 0.0 || l < fit_info.mod_depth_bounds_.first || l > fit_info.mod_depth_bounds_.second ) {
		// Object containing the best parameters so far
		struct {
			Real sum_of_squares_ = 0.0;
			Real d_ = 0.0;
			Real k_ = 0.0;
			Real l_ = 0.0;
		} best_stats;
		// Only look for d if we are also optimizing the dimensionality (e.g. for membrane proteins)
		Real min_d = ( optimize_dim ) ? 1.5 : 3.0;
		Real max_d = ( optimize_dim ) ? 3.5 : 3.0;
		Real const & min_l = fit_info.mod_depth_bounds_.first;
		Real const & max_l = fit_info.mod_depth_bounds_.second;
		// Outer loop: optimize dimensionality
		for ( Real xd = min_d; xd <= max_d; xd += 0.1 ) {
			// Inner loop" optimize modulation depth
			for ( Real xl = min_l; xl <= max_l; xl += 0.01 ) {
				Real xk = fit_k( exp_trace, trace, fit_info, xd, xl );
				auto temp_trace = fit_trace( exp_trace, trace, xd, xk, xl );
				Real trace_score = sum_of_squares( exp_trace, temp_trace );
				// Keep these parameters if they are better than anything so far
				// Note that NaN can happen when fitting non-3D data, it isn't a bug, it's a mathematical reality.
				if ( std::isnan( std::abs( trace_score ) ) ) {
					continue;
				} else if ( trace_score <= best_stats.sum_of_squares_ || best_stats.sum_of_squares_ == 0.0 ) {
					best_stats.sum_of_squares_ = trace_score;
					best_stats.d_ = xd;
					best_stats.k_ = xk;
					best_stats.l_ = xl;
				}
			}
		}
		k = best_stats.k_;
		d = best_stats.d_;
		l = best_stats.l_;
	}
	// Optimize these parameters
	optimize_mod_depth( exp_trace, trace, fit_info, d, k, l, optimize_dim );
	//TR << "k: " << k << std::endl;
	//TR << "l: " << l << std::endl;
	//TR << "d: " << d << std::endl;
	fit_info.last_slope_ = k;
	fit_info.last_mod_depth_ = l;
	fit_info.last_dim_ = d;
	return fit_trace( exp_trace, trace, d, k, l );
}

/// @brief Convert an intramolecular trace to an intermolecular trace given background parameters
std::map< Real, Real >
Simulated4PDEERTrace::fit_trace(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & trace,
	Real const & d,
	Real const & k,
	Real const & l
) {
	std::map< Real, Real > output;
	for ( auto const & time_trace : exp_trace ) {
		Real const & time = time_trace.first;
		output[ time ] = std::exp( k * pow( std::abs( time ), d / 3.0 ) ) * ( 1.0 - l * ( 1.0 - trace.at( time ) ) );
	}
	return output;
}

/// @brief Determine the residuals between an experimental trace and a simulated trace
Real
Simulated4PDEERTrace::sum_of_squares(
	std::map< Real, Real > const & exp_trace,
	Real const & noise // = 1.0
) {
	return sum_of_squares( exp_trace, trace_, noise );
}

/// @brief Determine the residuals between an experimental trace and a simulated trace
Real
Simulated4PDEERTrace::sum_of_squares(
	std::map< Real, Real > const & exp_trace,
	std::map< Real, Real > const & sim_trace,
	Real const & noise // = 1.0
) {
	Real output = 0.0;
	debug_assert( exp_trace.size() == sim_trace.size() );
	for ( auto const & time_trace : exp_trace ) {
		Real const & exp_time = time_trace.first;
		Real const & exp_signal = time_trace.second;
		Real sim_signal = 0.0;
		// Herein lies the danger of indexing a map using a Real: sometimes it doesn't exist
		if ( sim_trace.find( exp_time ) == sim_trace.end() ) {
			Real closest_proximity = 0.0;
			for ( auto const & sim_time_trace : sim_trace ) {
				Real proximity = std::abs( exp_time - sim_time_trace.first );
				if ( closest_proximity == 0.0 || proximity < closest_proximity ) {
					closest_proximity = proximity;
					sim_signal = sim_time_trace.second;
				}
			}
		} else {
			sim_signal = sim_trace.at( exp_time );
		}
		output += pow( ( exp_signal - sim_signal ) / noise, 2 );
	}
	// TR << "Sum of squares: " << output / trace.size() << std::endl;
	return output / sim_trace.size();
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
