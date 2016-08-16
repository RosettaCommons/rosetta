// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/SmoothPolicy.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/SmoothPolicy.hh>

// C/C++ headers
#include <cmath>

// Utility headers
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/minmax.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>

//Auto Headers
namespace protocols {
namespace nonlocal {

typedef utility::vector1<core::Real> ScoreList;

/// @brief Convenience method for retrieving the single highest-weighted sample
core::Size make_selection(numeric::random::WeightedReservoirSampler<core::Size>* sampler) {
	assert(sampler);
	utility::vector1<core::Size> results;
	sampler->samples(&results);
	return results[1];
}

SmoothPolicy::SmoothPolicy(core::fragment::FragSetCOP fragments)
: Policy(fragments) {}

core::Size SmoothPolicy::choose(const core::fragment::Frame& frame,
	const core::pose::Pose& pose) {
	using core::Size;
	using numeric::random::WeightedReservoirSampler;
	assert(frame.nr_frags() > 0);

	ScoreList scores;
	scorer_.score(frame, pose, scores);

	WeightedReservoirSampler<Size> sampler(1);
	for ( Size i = 1; i <= scores.size(); ++i ) {
		double score = scores[i];
		double fitness = std::sqrt(scorer_.cutoff() - score);

		if ( score < scorer_.cutoff() ) {
			sampler.consider_sample(i, fitness);
		}
	}

	// If no candidates met the score threshold, return the best scoring fragment.
	// Otherwise, randomly choose among the candidates by score.
	if ( sampler.num_considered() == 0 ) {
		return utility::argmin(scores);
	} else {
		return make_selection(&sampler);
	}
}

}  // namespace nonlocal
}  // namespace protocols
