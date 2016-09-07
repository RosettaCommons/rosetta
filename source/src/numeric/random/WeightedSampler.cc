// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/WeightedSampler.cc
///
/// @brief  a class to generate random positive integers using a given set of weights
/// @author Colin A. Smith


// Unit header or inline function header
#include <numeric/random/WeightedSampler.hh>

// Other numeric headers or inline function headers
#include <numeric/random/random.hh>

// External library headers
#include <utility/exit.hh>

// C++ headers

// Operating system headers

// Forward declarations


namespace numeric {
namespace random {

// All of the following should be ordered to match the header

// Static (class) constants

// Static (class) variables

// Private _inline_ methods -- and definitions

// Methods


WeightedSampler::WeightedSampler() :
	cumulative_distribution_valid_(false)
{
}

WeightedSampler::WeightedSampler(
	numeric::Size num_weights
) :
	weights_(num_weights, 0),
	cumulative_distribution_valid_(false)
{
}

WeightedSampler::WeightedSampler(
	utility::vector1<numeric::Real> const & weights
) :
	weights_(weights),
	cumulative_distribution_valid_(false)
{
}

WeightedSampler::WeightedSampler(
	WeightedSampler const & weighted_sampler
)
{
	*this = weighted_sampler;
}


WeightedSampler::~WeightedSampler() = default;


WeightedSampler &
WeightedSampler::operator=( WeightedSampler const & ) = default;


// void owning_ptr_acquire(WeightedSampler * p)
// {
// }


// void owning_ptr_release(WeightedSampler * p)
// {
// }

numeric::Size
WeightedSampler::random_sample(
	numeric::Real randnum
) const {
	assert(randnum >= 0);
	assert(randnum <= 1);

	if  ( !cumulative_distribution_valid_ ) update_cumulative_distribution();

	for ( numeric::Size i = 1; i <= cumulative_distribution_.size(); ++i ) {
		if ( cumulative_distribution_[i] && cumulative_distribution_[i] >= randnum ) return i;
	}

	return cumulative_distribution_.size();
}

numeric::Size
WeightedSampler::random_sample(
	numeric::random::RandomGenerator & rg
) const {
	return random_sample(rg.uniform());
}

void
WeightedSampler::update_cumulative_distribution() const {
	runtime_assert(weights_.size());

	cumulative_distribution_.resize(weights_.size());

	numeric::Real weight_sum(0);

	for ( numeric::Size i = 1; i <= weights_.size(); ++i ) {
		assert(weights_[i] >= 0);
		weight_sum += weights_[i];
	}

	cumulative_distribution_[1] = weights_[1]/weight_sum;

	for ( numeric::Size i = 2; i < weights_.size(); ++i ) {
		cumulative_distribution_[i] = cumulative_distribution_[i-1] + weights_[i]/weight_sum;
	}

	cumulative_distribution_.back() = 1;

	cumulative_distribution_valid_ = true;
}


} // namespace random
} // namespace numeric
