// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERMiscMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERMiscMethod_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERMiscMethod_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

class DEERMiscMethod : public DEERDistanceDistribution {
public:

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) override;

	/// @brief Set mode being used to score
	/// @param  val: Mode with which to score
	void
	mode(
		std::string const & val
	);

	/// @brief Return the mode being used to score
	/// @return Mode being used to score
	std::string
	mode() const;

private:

	/// @brief  Scoring modes available
	/// @detail (error thrown if mode assigned isn't among these)
	utility::vector1< std::string > modes_ = { "CONDITIONAL",
"BHATTACHARYYA", "HELLINGER" };

	/// @brief Scoring mode
	std::string mode_ = "CONDITIONAL";

};

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
