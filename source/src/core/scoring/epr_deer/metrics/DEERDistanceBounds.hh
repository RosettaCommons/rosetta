// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/DEERData.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  These classes contain the base and derived types for various DEER data containers.
///      The DEERData parent class stores generic information. The DEERDistanceBounds type
///      stores a distance value of interest and evaluates as a harmonic function. The
///      DEERDistanceDistribution type store the data as a probability distribution and
///      tries to maximize overlap. The DEERDecayData type stores the raw data and tries
///      to recapitulate it from the simulated distribution
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceBounds_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceBounds_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.fwd.hh>
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

/// @brief  Derived class for storing the data as a bounded function.
/// @detail Contains the upper and lower bounds, as well as steepness.
///     If the average simulated distance falls within these bounds the score is zero
///     If the average falls outside, it is evaluated via the steepness
///     For example: ( ( lb - d ) / s ) ^2 OR ( ( d  - ub ) / s ) ^2
///     lb: Lower bound
///     ub: Upper bound
///     d: Avg distance
///     s: Steepness
class DEERDistanceBounds : public DEERData {
public:

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) override;

	/// @brief  Returns the lower and upper distance bounds
	/// @return Pair of lower and upper bounds
	std::pair< Real, Real >
	bounds() const;

	/// @brief  Returns the step / steepness of the scoring function
	/// @return Steepness value
	Real const &
	step() const;

	/// @brief Sets the lower and upper bounds
	/// @param lo: Lower bound
	/// @param hi: Upper bound
	void
	bounds(
		Real const & lo,
		Real const & hi
	);

	/// @brief Sets the step / steepness
	/// @param step: Step
	void
	step(
		Real const & step
	);

private:

	/// @brief Lower bounds
	Real lo_ = 0.0;

	/// @brief Upper bound
	Real hi_ = 0.0;

	/// @brief Step
	Real step_ = 1.0;

};
} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
