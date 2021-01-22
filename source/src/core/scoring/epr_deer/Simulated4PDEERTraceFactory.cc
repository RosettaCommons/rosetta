// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/Simulated4PDEERTraceFactory.cc
/// @brief  Class that simulates the observable signal resulting from electron dipolar coupling
/// @details The dipolar coupling signal between two electrons can be simulated from a distance
///       distribution. This class does that simulation by pulling values from the datacache
///       container objects (DEERData) and effectively storing a vector with values of interest.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>
#include <core/scoring/epr_deer/util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

#include <numeric/constants.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix_operations.hh>
#include <numeric/MathVector.hh>
#include <numeric/nls/lmmin.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <map>
#include <tuple>
#include <cmath>
#include <cstdlib>

namespace core {
namespace scoring {
namespace epr_deer {

// Necessary for LM method below
using TriVec = std::tuple< utility::vector1< Real >,
	utility::vector1< Real >, utility::vector1< Real > >;

using TriVecC = std::tuple< utility::vector1< Real > const,
	utility::vector1< Real > const, utility::vector1< Real > const >;

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.Simulated4PDEERTraceFactory" );

/// @brief Function for adding non-3D background coupling to DEER trace
/// @param  sim_trace: Simulated DEER trace (intramolecular)
/// @param  time_pts: Time points for DEER trace
/// @param  depth: Modulation deoth
/// @param  slope: Background slope
/// @param  dim: Background dimensionality
/// @detail Note that this is a "friend" function, not part of the factory
/// @detail Also note that all values need to pass through sigmoid function
utility::vector1< Real >
add_bckg(
	utility::vector1< Real > const & sim_trace,
	utility::vector1< Real > const & time_pts,
	Real const & depth,
	Real const & slope,
	Real const & dim
) {
	Real const true_depth = sigmoid( depth, mindepth_, maxdepth_ );
	Real const true_slope = sigmoid( slope, minslope_, maxslope_ );
	Real const true_dim = sigmoid( dim, mindim_, maxdim_ );

	utility::vector1< Real > output( sim_trace.size(), 0.0 );
	for ( Size i = 1; i <= sim_trace.size(); ++i ) {
		Real const bckg = std::exp( -1 * pow( abs( time_pts[ i ] ) * true_slope,
			true_dim / 3.0 ) );
		output[ i ] = bckg * ( 1 - true_depth * ( 1 - sim_trace[ i ] ) );
	}
	return output;
}

/// @brief Unpacks data passed to lmmin function
/// @param  data: The data passed to lmmin
/// @param  time_points: Time point vector to assign
/// @param  exp_trace: Experimental DEER trace to assign
/// @param  sim_trace: Simulated DEER trace to assign
void
unpack(
	void const *data,
	utility::vector1< Real > & time_pts,
	utility::vector1< Real > & exp_trace,
	utility::vector1< Real > & sim_trace
) {
	std::tie( time_pts, exp_trace, sim_trace ) =
		*( static_cast< TriVec * >( const_cast< void* >( data ) ) );
}

/// @brief LM function for optimizing non-3D background coupling
/// @param  par: C-style parameter array
/// @param  m_dat: Number of data points
/// @param  data: Generic data pointer object that needs to be unpacked
/// @param  fvec: Vector of residuals used by lmmin()
/// @param  null (not used)
/// @detail Note that this is a "friend" function, not part of the factory
void
minimize_sumsq(
	Real const *par,
	int m_dat,
	void const *data,
	Real *fvec,
	int * // info
) {

	// First unpack our vectors
	utility::vector1< Real > time_pts, exp_trace, sim_trace;
	unpack( data, time_pts, exp_trace, sim_trace );

	// Then simulate the data and take the sum of squared residuals
	auto const sim_bckg = add_bckg( sim_trace, time_pts, par[ 0 ], par[ 1 ],
		par[ 2 ] );
	for ( int i = 0; i < m_dat; ++i ) {
		fvec[ i ] = pow( sim_bckg[ i + 1 ] - exp_trace[ i + 1 ], 2 );
	}
}

/// @brief LM function for optimizing 3D background coupling
/// @param  par: C-style parameter array
/// @param  m_dat: Number of data points
/// @param  data: Generic data pointer object that needs to be unpacked
/// @param  fvec: Vector of residuals used by lmmin()
/// @param  null (not used)
/// @detail Note that this is a "friend" function, not part of the factory
void
minimize_sumsq_nodim(
	Real const *par,
	int m_dat,
	void const *data,
	Real *fvec,
	int * // info
) {

	// First unpack our vectors
	utility::vector1< Real > time_pts, exp_trace, sim_trace;
	unpack( data, time_pts, exp_trace, sim_trace );

	// Then simulate the data and take the sum of squared residuals
	auto const sim_bckg = add_bckg( sim_trace, time_pts, par[ 0 ], par[ 1 ],
		dim3d_ );
	for ( int i = 0; i < m_dat; ++i ) {
		fvec[ i ] = pow( sim_bckg[ i + 1 ] - exp_trace[ i + 1 ], 2 );
	}
}

/// @brief  More comprehensive constructor that gets everything in order
/// @param  exp_trace: Experimental DEER trace for comparison
/// @param  time_pts: Time points for the experimental DEER trace
/// @param  bckg_type: Type of background. Can be "NONE", "3D", "NON-3D"
/// @param  max_dist: Maximum distance in angstroms
/// @param  n_bins: For discretization of angle-dependent dipolar coupling
Simulated4PDEERTraceFactory::Simulated4PDEERTraceFactory(
	utility::vector1< Real > const & exp_trace,
	utility::vector1< Real > const & time_pts,
	std::string const & bckg_type,
	Size const & bins_per_a,
	Size const & max_dist, // = 100
	Size const & n_bins // = 200
) {

	// Ensure that the two vectors are the same size
	if ( exp_trace.size() != time_pts.size() ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Sizes of time point vector and DEER trace do not match!" );
	}

	// Set the values as provided
	exp_data_ = exp_trace;
	time_pts_ = time_pts;
	bckg_type_ = bckg_type;
	bins_per_a_ = bins_per_a;
	max_dist_ = max_dist;
	n_bins_ = n_bins;

	// Calculate the kernel matrix used to convert distances to decays
	kernel_ = initialize_kernel( time_pts, bins_per_a, max_dist, n_bins );
}

/// @brief  Kernel matrix initialization method
/// @param  time_pts: Time points to be used in the kernel
/// @param  bins_per_a: Bins per angstrom
/// @param  max_dist: Maximum distance in Angstroms to be stored here
/// @param  n_bins: For discretization of angle-dependent dipolar coupling
/// @return Kernel matrix
/// @detail Mathematical details can be found in Ibanez and Jeschke, 2019
/// @detail   JMR (doi: 10.1016/j.jmr.2019.01.008)
numeric::MathMatrix< Real >
Simulated4PDEERTraceFactory::initialize_kernel(
	utility::vector1< Real > const & time_pts,
	Size const & bins_per_a,
	Size const & max_dist, // = 100
	Size const & n_bins // = 200
) const {

	// Initialize matrix with correct size
	Size const & n_rows = time_pts.size();
	Size const n_cols = bins_per_a * max_dist;
	numeric::MathMatrix< Real > kernel( n_rows, n_cols );

	// Iterate through rows (time) and columns (distance)
	for ( Size row = 0; row < n_rows; ++row ) {
		Real t = time_pts[ row + 1 ];
		for ( Size col = 0; col < n_cols; ++col ) {

			// Dipolar coupling constant
			Real const c_top = t * 326.98;
			Real const c_bot = pow( Real( bins_per_a * col + 1 ) / 10, 3 );
			Real const c = c_top / c_bot;
			Real val = 0.0;

			// Now integrate from zero to pi_half (angle-dependent coupling)
			for ( Size bin = 0; bin <= n_bins; ++bin ) {
				Real const a = bin * numeric::constants::d::pi_over_2 / n_bins;
				Real const coupling = 1 - 3 * pow( std::cos( a ), 2 );
				val += std::sin( a ) * std::cos( coupling ) * ( c );
			}

			// Assign result to position of interest in kernel
			kernel( row, col ) = val / ( n_bins + 1 );
		}
	}
	return kernel;
}

/// @brief  Process for converting distance distribution to DEER trace
/// @param  distr: DEER distribution stored as a map
/// @return DEER trace with all the bells and whistles attached
Simulated4PDEERTrace
Simulated4PDEERTraceFactory::trace_from_distr(
	std::map< Size, Real > const & distr
) {

	// Obtain the dipolar coupling decay data within a molecule
	auto const deer_trace_intra = kernel_mult( distr );

	// Calculate the intermolecular contribution by least-squares
	return Simulated4PDEERTrace( opt_bckg( deer_trace_intra ),
		deer_trace_intra, time_pts_, depth_, slope_, dim_ );
}

/// @brief Actual nitty-gritty code for applying the kernel matrix
/// @param  distr: DEER distribution stored as a map
/// @return Vector with intramolecular component
/// @detail Note that throughout this calculation, the bins per angstrom
/// @detail  and time points are assumed to be correct! This is not
/// @detail  checked (maybe a debug assert can be used?)
utility::vector1< Real >
Simulated4PDEERTraceFactory::kernel_mult(
	std::map< Size, Real > const & distr
) const {

	// Convert map to the right object
	numeric::MathVector< Real > d_vec( kernel_.get_number_cols(), 0.0 );
	for ( Size i = 0; i < kernel_.get_number_cols(); ++i ) {
		if ( distr.find( i ) != distr.end() ) {
			d_vec( i ) = distr.at( i );
			//  TR << "\t" << i << "\t" << distr.at( i ) << std::endl;
			// } else {
			//  TR << "\t" << i << "\t" << 0 << std::endl;
		}
	}

	// The actual conversion
	numeric::MathVector< Real > const trace = kernel_ * d_vec;

	// Convert to vector1 object with normalized signal
	Real maxval = *std::max_element( trace.begin(), trace.end() );

	// Edge cases if something doesn't work
	if ( maxval == 0.0 ) {
		TR.Error << "Simulated DEER trace is full of zeros!" << std::endl;
		return utility::vector1< Real >( trace.size(), 1.0 );
	} else if ( std::isinf( std::abs( maxval ) ) ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Maximum value is infinite!" );
	} else if ( std::isnan( std::abs( maxval ) ) ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Maximum value is NaN!" );
	}

	// Otherwise, proceed as planned and return the output
	utility::vector1< Real > output;
	std::transform( trace.begin(), trace.end(), std::back_inserter( output ),
		[&]( Real const & val ){ return val / maxval; } );
	return output;
}

/// @brief  Optimize the intermolecular background by least squares
/// @param  deer_trace: The intramolecular deer trace
/// @return vector with final DEER trace
/// @detail The function can be run with and without dimensionality
/// @detail  calculation. This is reserved for cases where membrane
/// @detail  proteins are studied and are not uniformly distributed
/// @detail  in three-dimensional space.
utility::vector1< Real >
Simulated4PDEERTraceFactory::opt_bckg(
	utility::vector1< Real > const & deer_trace
) {

	// Skip this step if no background calculation is necessary
	if ( bckg_type_ == "NONE" ) {
		depth_ = 0.0;
		slope_ = 0.0;
		dim_ = 0.0;
		return deer_trace;
	} else {

		// The "data" being passed to LM optimization
		numeric::nls::lm_status_struct status;
		TriVec tuple_data = std::make_tuple( time_pts_, exp_data_, deer_trace );
		auto * data = &tuple_data;
		utility::vector1< Real > params{ depth_, slope_ };

		// Number of iterations
		auto const n_iter = ( optimized_ ) ? runs_ : init_runs_;

		// The function itself depends on if a non-3D background is used
		if ( bckg_type_ != "3D" ) {
			params.push_back( dim_ );
			for ( Size i = 1; i <= n_iter; ++i ) {
				numeric::nls::lmmin( params.size(), &params[ 1 ],
					time_pts_.size(), ( const void* )data, minimize_sumsq,
					&status, numeric::nls::lm_printout_std );
			}

			// If dimensionality does not need to be optimized (fixed at 3.0)
		} else {
			for ( Size i = 1; i <= n_iter; ++i ) {
				numeric::nls::lmmin( params.size(), &params[ 1 ],
					time_pts_.size(), ( const void* )data, minimize_sumsq_nodim,
					&status, numeric::nls::lm_printout_std );
			}
		}

		// Note that parameters have been tuned; save the background
		// parameters and return the DEER trace
		optimized_ = true;
		depth_ = params[ 1 ];
		slope_ = params[ 2 ];
		if ( bckg_type_ != "3D" ) {
			dim_ = params[ 3 ];
		}

		return add_bckg( deer_trace, time_pts_, depth_, slope_, dim_ );
	}
}

/// @brief  Adds a time point to the saved experimental DEER data
/// @param  trace_datum: The datapoint (Y-axis)
/// @param  time_point: The datapoint's time (X-axis)
/// @detail NOT RECOMMENDED. The entire kernel matrix is reinitialized!
/// @detail This is insanely time-consuming!
void
Simulated4PDEERTraceFactory::append_data_and_recalculate(
	Real const & trace_datum,
	Real const & time_point
) {
	exp_data_.push_back( trace_datum );
	time_pts_.push_back( time_point );
	kernel_ = initialize_kernel( time_pts_, bins_per_a_, max_dist_, 200 );
}

/// @brief Returns experimental data
/// @return DEER data
utility::vector1< Real >
Simulated4PDEERTraceFactory::trace() const {
	return exp_data_;
}

/// @brief  Returns time points for experimental data
/// @return Time points for experimental data
utility::vector1< Real >
Simulated4PDEERTraceFactory::time_pts() const {
	return time_pts_;
}

/// @brief  Returns background type
/// @return Background type
std::string
Simulated4PDEERTraceFactory::bckg_type() const {
	return bckg_type_;
}

/// @brief  Returns bins per angstrom in the kernel matrix
/// @return Bins per angstrom in the kernel matrix
Size
Simulated4PDEERTraceFactory::bins_per_a() const {
	return bins_per_a_;
}

/// @brief Sets experimental data to new DEER trace
/// @param  val: New data to set
void
Simulated4PDEERTraceFactory::trace(
	utility::vector1< Real > const & val
) {
	exp_data_ = val;
}

/// @brief Sets time points to new data and recalculates matrix
/// @param  val: New data to set
void
Simulated4PDEERTraceFactory::time_pts(
	utility::vector1< Real > const & val
) {
	time_pts_ = val;
	kernel_ = initialize_kernel( time_pts_, bins_per_a_, max_dist_, n_bins_ );
}

/// @brief Sets background coupling type
/// @param  val: New background coupling
void
Simulated4PDEERTraceFactory::bckg_type(
	std::string const & val
) {
	bckg_type_ = val;
}

/// @brief  Sets bins per angstrom of kernel matrix and recalculates matrix
/// @param  val: New value for bins per angstrom
void
Simulated4PDEERTraceFactory::bins_per_a(
	Size const & val
) {
	bins_per_a_ = val;
	kernel_ = initialize_kernel( time_pts_, bins_per_a_, max_dist_, n_bins_ );
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
