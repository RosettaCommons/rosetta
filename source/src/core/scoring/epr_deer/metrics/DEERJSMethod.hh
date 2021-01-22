// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/epr_deer/metrics/DEERJSMethod.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERJSMethod_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERJSMethod_hh

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
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

class DEERJSMethod : public DEERDistanceDistribution {
public:

	/// @brief  Virtual function to evaluate score given a distribution
	/// @param  sim_histr: Simulated DEER distribution
	/// @return Freshly computed score
	Real
	get_score(
		std::map< Size, Real > const & sim_histr
	) override;

};

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
