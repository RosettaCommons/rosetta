// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh
/// @brief   Container for DEER experimental data and dataype-specific scoring function
/// @detail  This class stores the data as a probability distribution and
///        tries to maximize some function (see derived types).
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceDistribution_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERDistanceDistribution_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>

#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <string>
#include <map>
#include <tuple>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

/// Derived class that stores the entire distance distribution
/// Score is evaluated using the cross-entropy of the simulated from the experimental
class DEERDistanceDistribution : public DEERData {

public:

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	/// @detail Cross-entropy corresponds to the negative log-likelihood that the
	///     experimental distribution could have given rise to the simulated. This
	///     allows boltzmann weighting and/or Bayesian statistical inference from score
	///     Note that although confidence bands can be received as input, they are
	///     not currently used for this purpose. If you know an information-theoretic
	///     approach to using them, please contact me. I would live to incorporate
	///     that information here.
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) override;

	/// @brief Expand distance distribution into multiple distributions
	/// @param distr: Simulated distribution to expand
	/// @return Vector of distributions
	/// @detail Used for Bayesian analysis when backbone dynamics are not known
	///    and several possibilities must be taken into account
	utility::vector1< std::map< Size, Real > >
	expand_hists(
		std::map< Size, Real > const & distr
	);

	/// @brief Get range of possible P(r) values for experimental distribution
	/// @param  bin: Value of r
	/// @return Vector of P(r) for Bayesian scoring with confidence bands
	utility::vector1< Real >
	get_prs(
		Size const & bin
	) const;

	// SETTTERS

	/// @brief Sets the lower bound/confidence band for the distance distribution
	/// @param val: Lower bound or confidence band. Note: Can be negative!
	void
	lower_bound(
		std::map< Size, Real > const & val
	);

	/// @brief Sets the line of best fit
	/// @param val: Best fit of DEER distribution
	void
	best_fit(
		std::map< Size, Real > const & val
	);

	/// @brief Sets the upper bound/confidence band for the distance distribution
	/// @param val: Lower bound or confidence band.
	void
	upper_bound(
		std::map< Size, Real > const & val
	);

	/// @brief Set whether confidence bands will be used when calculating score
	/// @param  val: Whether to use confidence bands
	void
	bounds(
		bool const & val
	);

	/// @brief Set whether to compute the reverse metric
	/// @param  val: Whether to compute the reverse metric
	void
	reverse(
		bool const & val
	);

	/// @brief Set whether to use confidence bands
	/// @param  val: Whether to use confidence bands
	void
	bb(
		bool const & val
	);

	/// @brief Set whether to use a single distance
	/// @param  val: Whether to use a single distance
	void
	singleval(
		bool const & val
	);

	/// @brief Set whether to compute the integral
	/// @param  val: Whether to compute the integral
	void
	integral(
		bool const & val
	);

	// GETTERS

	/// @brief  Returns the lower bound/confidence band for the distribution
	/// @return The lower bound/confidence band for the distribution
	std::map< Size, Real > const &
	lower_bound() const;

	/// @brief Returns the best fit of the distance distribution
	/// @brief The best fit of the distance distribution
	std::map< Size, Real > const &
	best_fit() const;

	/// @brief  Returns the upper bound/confidence band for the distribution
	/// @return The upper bound/confidence band for the distribution
	std::map< Size, Real > const &
	upper_bound() const;

	/// @brief  Returns whether confidence bands are being used
	/// @return Whether confidence bands are being used
	bool
	bounds() const;

	/// @brief  Returns whether the reverse metric is being used
	/// @return Whether the reverse metric is being used
	bool
	reverse() const;

	/// @brief  Returns whether backbone expansion is being used
	/// @return Whether backbone expansion is being used
	bool
	bb() const;

	/// @brief  Returns whether calculation of a single distance is being used
	/// @return Whether calculation of a single distance is being used
	bool
	singleval() const;

	/// @brief  Returns whether calculation of the integral is being used
	/// @return Whether calculation of the integral is being used
	bool
	integral() const;

protected:

	/// @brief Best fit / distance distribution
	std::map< Size, Real > distr_ = {};

	/// @brief Lower bound for distribution (95% confidence interval)
	std::map< Size, Real > lower_bound_ = {};

	/// @brief Upper bound for distribution (95% confidence interval)
	std::map< Size, Real > upper_bound_ = {};

	/// @brief Whether the calculation method uses confidence bands
	bool bounds_ = false;

	/// @brief Whether the reverse metric is computed
	bool rev_ = false;

	/// @brief Whether backbone expansion is done (for Bayesian scoring)
	bool bb_ = false;

	/// @brief Whether a single value is being used for scoring
	bool singleval_ = false;

	/// @brief Whether the integral is computed
	bool integral_ = false;

private:

	/// @brief  The standard deviation of the experimental distribution
	Real exp_stdev_ = 0.0;
};

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
