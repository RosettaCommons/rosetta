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

#include <core/scoring/epr_deer/Simulated4PDEERTrace.fwd.hh>

#include <core/scoring/epr_deer/DEERData.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

#include <map>

namespace core {
namespace scoring {
namespace epr_deer {

class Simulated4PDEERTrace {

public:

	/// @brief Default constructor
	Simulated4PDEERTrace();

	/// @brief Custom constructor
	Simulated4PDEERTrace(
		std::map< Real, Real > const & exp_trace,
		FittingInfo & fit_info,
		std::map< Size, Real > const & sim_histogram,
		Size const & bins_per_a
	);

	/// @brief Bracket operator for accessing the DEER trace
	Real
	operator[](
		Real const & t
	) const;

	/// @brief Main function for simulating 4-pulse DEER trace
	/// @detail This function first simulates the intramolecular dipolar coupling
	///     form factor using the evolution kernel described in
	///     "DEER Distance Measurements on Proteins" by Jeschke (2012).
	///     It then finds the best background fit for the experimental data,
	///     if necessary (3D by default, but can be set to "NONE" or "NON-3D")
	///     Background-fitted DEER trace is saved and returned.
	std::map< Real, Real >
	simulate_decay(
		std::map< Real, Real > const & exp_trace,
		FittingInfo & fit_info,
		std::map< Size, Real > const & sim_histogram,
		Size const & bins_per_a
	);

	/// @brief Deterine the intramolecular DEER trace
	std::map< Real, Real >
	calc_4pdeer(
		std::map< Real, Real > const & exp_trace,
		FittingInfo const & fit_info,
		std::map< Size, Real > const & sim_histogram,
		Size const & bins_per_a
	);

	/// @brief Fit the background decay using linear regression
	Real
	fit_k(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & trace,
		FittingInfo & fit_info,
		Real const & d,
		Real const & l
	);

	/// @brief Optimize the background function of the DEER trace: general function
	std::map< Real, Real >
	optimize_bckg(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & trace_no_bckg,
		FittingInfo & fit_info
	);

	/// @brief Optimize the modulation depth using simple gradient descent
	void
	optimize_mod_depth(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & trace,
		FittingInfo & fit_info,
		Real & d,
		Real & k,
		Real & l,
		bool fit_dim // = false
	);

	/// @brief Optimize the background function of the DEER trace: detailed function
	std::map< Real, Real >
	optimize_bckg(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & trace,
		FittingInfo & fit_info,
		bool const & optimize_dim
	);

	/// @brief Convert an intramolecular trace to an intermolecular trace given background parameters
	std::map< Real, Real >
	fit_trace(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & trace,
		Real const & d,
		Real const & k,
		Real const & l
	);

	/// @brief Determine the residuals between an experimental trace and a simulated trace
	Real
	sum_of_squares(
		std::map< Real, Real > const & exp_trace,
		Real const & noise = 1.0
	);

	/// @brief Determine the residuals between an experimental trace and a simulated trace
	Real
	sum_of_squares(
		std::map< Real, Real > const & exp_trace,
		std::map< Real, Real > const & sim_trace,
		Real const & noise = 1.0
	);

private:

	std::map< Real, Real > trace_;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core
