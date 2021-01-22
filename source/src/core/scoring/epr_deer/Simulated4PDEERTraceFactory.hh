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

#ifndef INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_hh
#define INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_hh

#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.fwd.hh>

#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>
#include <core/types.hh>

#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>

#include <utility/vector1.hh>

#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {

//////////////////////////////////////////////////
// OBJECTS AND FUNCTIONS RELATED TO OPTIMIZATION

/// @brief Minimum modulation depth (passed to sigmoid function)
static Real const mindepth_ = 0.01;

/// @brief Maximum modulation depth (passed to sigmoid function)
static Real const maxdepth_ = 0.9;

/// @brief Minimum background slope (passed to sigmoid function)
static Real const minslope_ = 0.0;

/// @brief Maximum background slope (passed to sigmoid function)
static Real const maxslope_ = 1.0;

/// @brief Minimum dimensionality (passed to sigmoid function)
static Real const mindim_ = 2.0;

/// @brief Maximum dimensionality (passed to sigmoid function)
static Real const maxdim_ = 3.5;

/// @brief Default dimensionality to yield 3.0 (passed to sigmoid function)
static Real const dim3d_ = 0.69314717;

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
);

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
);

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
);

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
);

class Simulated4PDEERTraceFactory {

public:

	/// @brief  Default constructor
	Simulated4PDEERTraceFactory() = default;

	/// @brief  More comprehensive constructor that gets everything in order
	/// @param  exp_trace: Experimental DEER trace for comparison
	/// @param  time_pts: Time points for the experimental DEER trace
	/// @param  bckg_type: Type of background. Can be "NONE", "3D", "NON-3D"
	/// @param  max_dist: Maximum distance in angstroms
	/// @param  n_bins: For discretization of angle-dependent dipolar coupling
	Simulated4PDEERTraceFactory(
		utility::vector1< Real > const & exp_trace,
		utility::vector1< Real > const & time_pts,
		std::string const & bckg_type,
		Size const & bins_per_a,
		Size const & max_dist = 100,
		Size const & n_bins = 200
	);

	/// @brief  Kernel matrix initialization method
	/// @param  time_pts: Time points to be used in the kernel
	/// @param  bins_per_a: Bins per angstrom
	/// @param  max_dist: Maximum distance in Angstroms to be stored here
	/// @param  n_bins: For discretization of angle-dependent dipolar coupling
	/// @return Kernel matrix
	/// @detail Mathematical details can be found in Ibanez and Jeschke, 2019
	/// @detail   JMR (doi: 10.1016/j.jmr.2019.01.008)
	numeric::MathMatrix< Real >
	initialize_kernel(
		utility::vector1< Real > const & time_pts,
		Size const & bins_per_a,
		Size const & max_dist = 100,
		Size const & n_bins = 200
	) const;

	/// @brief  Process for converting distance distribution to DEER trace
	/// @param  distr: DEER distribution stored as a map
	/// @return DEER trace with all the bells and whistles attached
	Simulated4PDEERTrace
	trace_from_distr(
		std::map< Size, Real > const & distr
	);

	/// @brief Actual nitty-gritty code for applying the kernel matrix
	/// @param  distr: DEER distribution stored as a map
	/// @return Vector with intramolecular component
	/// @detail Note that throughout this calculation, the bins per angstrom
	/// @detail  and time points are assumed to be correct! This is not
	/// @detail  checked (maybe a debug assert can be used?)
	utility::vector1< Real >
	kernel_mult(
		std::map< Size, Real > const & distr
	) const;

	/// @brief  Optimize the intermolecular background by least squares
	/// @param  deer_trace: The intramolecular deer trace
	/// @return vector with final DEER trace
	/// @detail The function can be run with and without dimensionality
	/// @detail  calculation. This is reserved for cases where membrane
	/// @detail  proteins are studied and are not uniformly distributed
	/// @detail  in three-dimensional space.
	utility::vector1< Real >
	opt_bckg(
		utility::vector1< Real > const & deer_trace
	);

	/// @brief  Adds a time point to the saved experimental DEER data
	/// @param  trace_datum: The datapoint (Y-axis)
	/// @param  time_point: The datapoint's time (X-axis)
	/// @detail NOT RECOMMENDED. The entire kernel matrix is reinitialized!
	/// @detail This is insanely time-consuming!
	void
	append_data_and_recalculate(
		Real const & trace_datum,
		Real const & time_point
	);

	/// @brief Returns experimental data
	/// @return DEER data
	utility::vector1< Real >
	trace() const;

	/// @brief  Returns time points for experimental data
	/// @return Time points for experimental data
	utility::vector1< Real >
	time_pts() const;

	/// @brief  Returns background type
	/// @return Background type
	std::string
	bckg_type() const;

	/// @brief  Returns bins per angstrom in the kernel matrix
	/// @return Bins per angstrom in the kernel matrix
	Size
	bins_per_a() const;

	/// @brief Sets experimental data to new DEER trace
	/// @param  val: New data to set
	void
	trace(
		utility::vector1< Real > const & val
	);

	/// @brief Sets time points to new data
	/// @param  val: New data to set
	void
	time_pts(
		utility::vector1< Real > const & val
	);

	/// @brief Sets background coupling type
	/// @param  val: New background coupling
	void
	bckg_type(
		std::string const & val
	);

	/// @brief  Sets bins per angstrom of kernel matrix and recalculates matrix
	/// @param  val: New value for bins per angstrom
	void
	bins_per_a(
		Size const & val
	);

private:

	/// @brief  Experimental data
	utility::vector1< Real > exp_data_;

	/// @brief  Time points
	utility::vector1< Real > time_pts_;

	/// @brief  Background type ("3D", "NON-3D", "NONE")
	std::string bckg_type_;

	/// @brief  Bins per angstrom in kernel matrix
	Size bins_per_a_;

	Size max_dist_;

	/// @brief  Number of bins for magnetic field vector angle calculation
	Size n_bins_ = 200;

	/// @brief Matrix that converts distribution to intramolecular DEER trace
	numeric::MathMatrix< Real > kernel_;

	/// @brief  Last saved modulation depth
	Real depth_ = 0.25;

	/// @brief  Last saved background slope
	Real slope_ = 1e-4;

	/// @brief  Last saved dimensionality
	/// @detail Default value chosen because sigmoid( 0.6931, 2, 3.5 ) = 3
	/// @detail See sigmoid function in util.cc for details
	Real dim_ = dim3d_;

	/// @brief  To know if more LM iterations need to be run
	bool optimized_ = false;

	/// @brief How many times to run LM the first time
	Size init_runs_ = 10;

	/// @brief How many times to run LM subsequent times
	Size runs_ = 2;

};

/*

/// @brief Struct for fitting DEER traces. Copies are stored in each DEERDecayData object
struct FittingInfo {
Real last_slope_ = 0.0;
Real last_mod_depth_ = 0.0;
Real last_dim_ = 3.0;
Real time_pts_sqd_ = 0.0;
std::map< Size, std::map< Size, Real > > spin_map_; // distance outer, time inner
std::pair< Real, Real > mod_depth_bounds_ = std::make_pair( 0.02, 0.75 );
std::string bckg_ = "3D";
};

class Simulated4PDEERTraceFactory {

public:

/// @brief Default constructor
Simulated4PDEERTraceFactory();

/// @brief Custom constructor
Simulated4PDEERTraceFactory(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
metrics::FittingInfo & fit_info,
std::map< Size, Real > const & sim_histogram,
Size const & bins_per_a
);

/// @brief Bracket operator for accessing the DEER trace
Size
operator[](
Size const & t
) const;

/// @brief Main function for simulating 4-pulse DEER trace
/// @detail This function first simulates the intramolecular dipolar coupling
///     form factor using the evolution kernel described in
///     "DEER Distance Measurements on Proteins" by Jeschke (2012).
///     It then finds the best background fit for the experimental data,
///     if necessary (3D by default, but can be set to "NONE" or "NON-3D")
///     Background-fitted DEER trace is saved and returned.
utility::vector1< Real >
simulate_decay(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
metrics::FittingInfo & fit_info,
std::map< Size, Real > const & sim_histogram,
Size const & bins_per_a
);

/// @brief Deterine the intramolecular DEER trace
utility::vector1< Real >
calc_4pdeer(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
metrics::FittingInfo const & fit_info,
std::map< Size, Real > const & sim_histogram,
Size const & bins_per_a
);

/// @brief Fit the background decay using linear regression
Real
fit_k(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo & fit_info,
Real const & d,
Real const & l
);

/// @brief Optimize the background function of the DEER trace: general function
utility::vector1< Real >
optimize_bckg(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo & fit_info
);

/// @brief Optimize the modulation depth using simple gradient descent
void
optimize_mod_depth(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo & fit_info,
Real & d,
Real & k,
Real & l,
bool fit_dim // = false
);

/// @brief Optimize the background function of the DEER trace: detailed function
utility::vector1< Real >
optimize_bckg(
utility::vector1< Real > const & exp_trace,
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo & fit_info,
bool const & optimize_dim
);

/// @brief Convert an intramolecular trace to an intermolecular trace given background parameters
utility::vector1< Real >
fit_trace(
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
Real const & d,
Real const & k,
Real const & l
);

utility::vector1< Real >
fit_trace(
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo const & info
);

utility::vector1< Real >
return_intra(
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
metrics::FittingInfo const & info
);

utility::vector1< Real >
return_intra(
std::map< Size, Real > const & time_pts,
utility::vector1< Real > const & trace,
Real const & d,
Real const & k,
Real const & l
);

/// @brief Determine the residuals between an experimental trace and a simulated trace
Real
sum_of_squares(
utility::vector1< Real > const & exp_trace,
Real const & noise = 1.0
);

/// @brief Determine the residuals between an experimental trace and a simulated trace
Real
sum_of_squares(
utility::vector1< Real > const & exp_trace,
utility::vector1< Real > const & sim_trace,
Real const & noise = 1.0
);

private:

utility::vector1< Real > trace_;
utility::vector1< Real > time_pts_;

Real last_slope_ = 0.0;
Real last_mod_depth_ = 0.0;
Real last_dim_ = 3.0;
Real time_pts_sqd_ = 0.0;
std::map< Size, std::map< Size, Real > > spin_map_; // distance outer, time inner
std::pair< Real, Real > mod_depth_bounds_ = std::make_pair( 0.02, 0.75 );
std::string bckg_ = "3D";



};

*/

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_hh
