// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERDecayData.hh
/// @brief   Container for DEER experimental decay data
/// @detail  This class stores the raw data and tries
///        to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERDecayData_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERDecayData_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.fwd.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.hh>
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

/// @brief Derived class that stores DEER decay data
class DEERDecayData : public DEERData {
public:

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	/// @detail A DEER trace is calculated from the simulated histogram.
	///  For details, please read del Alamo Biophysical Journal 2020.
	///  The score is the average sum-of-squares residuals between the
	///  two. Note that there is a slight noise-dependence to this;
	///  noisier data will always return higher scores.
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) override;

	/// @brief Calculate sum of squares of simulated DEER trace
	/// @param  trace1: Simulated DEER trace
	/// @param  trace2: Reference DEER trace; experimental by default
	/// @param  normalize: Normalize by number of time points; recommended
	/// @return Sum of squared residuals
	Real
	sum_of_squares(
		utility::vector1< Real > const & sim_trace,
		bool const & normalize = true
	) const;

	/// @brief Initialize DEER factory object
	/// @param trace: Experimental DEER trace
	/// @param tiem_pts: Time points corresponding to DEER trace
	void
	init_factory(
		utility::vector1< Real > const & trace,
		utility::vector1< Real > const & time_pts
	);

	/// @brief Returns DEER trace factory object
	/// @return Factory object
	Simulated4PDEERTraceFactory &
	factory();

	/// @brief  Returns the noise from the imaginary component
	/// @return Noise level
	Real const &
	noise() const;

	/// @brief  Returns whether standard deviation is a fitting parameter
	/// @return Boolean
	bool const &
	fit_stdev() const;

	/// @brief Returns background type
	/// @return Intermolecular coupling background type
	std::string
	bckg() const;

	/// @brief Returns maximum distance for kernel calculation
	/// @return Maximum distance for kernel calculation
	Size
	max_dist() const;

	/// @brief  Sets DEER trace factory object
	/// @param  factory: Factory object
	void
	factory(
		Simulated4PDEERTraceFactory const & val
	);

	/// @brief Sets the noise from the imaginary component
	/// @param val: Noise level
	void
	noise(
		Real const & val
	);

	/// @brief Sets if standard deviation is varied as a parameter
	/// @param val: Boolean
	void
	fit_stdev(
		bool const & val
	);

	/// @brief Sets background type
	/// @param  val: Type of intermolecular coupling background
	void
	bckg(
		std::string const & val
	);

	/// @brief Sets maximum distance for kernel calculation
	/// @param val: Maximum distance to set
	void
	max_dist(
		Size const & val
	);

private:

	/// @brief DEER trace factory object
	Simulated4PDEERTraceFactory factory_;

	/// @brief Noise level in experimental DEER data
	Real noise_ = 1.0;

	/// @brief Whether standard deviation should be varied as a parameter
	bool fit_stdev_ = false;

	/// @brief Intermolecular background coupling type
	std::string bckg_ = "3D";

	/// @brief Maximum distance considered for kernel calculation
	Size max_dist_ = 100;

};

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
