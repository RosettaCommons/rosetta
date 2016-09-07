// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/reservoir_sample.hh
/// @brief  Randomly select the best N elements from a stream of elements using one
/// pass over a dataset.
/// @details
/// @author James Thompson


#ifndef INCLUDED_numeric_random_reservoir_sample_hh
#define INCLUDED_numeric_random_reservoir_sample_hh

// Unit headers
#include <numeric/types.hh>
#include <numeric/random/random.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ Headers

namespace numeric {
namespace random {

/// @brief Returns the probability that the Nth value in a sequence
/// should be accepted using the reservoir sampling criterion.
/// @details If we've seen N values and we want to keep K of them,
/// the probability of the Nth value being accepted is min(K/N,1.0).
inline numeric::Real
reservoir_sample_accept_prob(
	numeric::Size n_wanted,
	numeric::Size n_seen
) {
	Real const accept_prob(
		static_cast< Real > (n_wanted) / static_cast< Real > (n_seen)
	);
	return std::min( accept_prob, 1.0 );
}

/// @brief Simple container for keeping K random values.
/// @details Values are stochastically preserved so that the probability of
/// any value existing is always 1 / N, where N is the number of values seen.
template< typename T >
class ReservoirSampler {
public:
	ReservoirSampler( numeric::Size const wanted ) :
		n_seen_(0),
		n_wanted_( wanted )
	{}

	~ReservoirSampler() = default;

	void add_value( T const & val ) {
		++n_seen_;
		if ( n_vals() < n_wanted() ) {
			values_.push_back( val );
		} else {
			Real const accept_prob(
				reservoir_sample_accept_prob( n_wanted(), n_seen() )
			);

			if ( numeric::random::uniform() <= accept_prob ) {
				Size const idx = numeric::random::random_range( 1, n_vals() );
				values_[ idx ] = val;
			}
		}
	} // add_value

	numeric::Size n_vals() const {
		return values_.size();
	}

	numeric::Size n_wanted() const {
		return n_wanted_;
	}

	numeric::Size n_seen() const {
		return n_seen_;
	}

	utility::vector1< T > const & values() const {
		return values_;
	}

private:

	ReservoirSampler(); // Must pass n_wanted to constructor

	utility::vector1< T > values_;
	numeric::Size n_seen_;
	numeric::Size const n_wanted_;
};

template< typename T >
utility::vector1< T >
reservoir_sample(
	utility::vector1< T > const & vec,
	numeric::Size n_wanted,
	RandomGenerator & rg = numeric::random::rg()
) {
	using numeric::Real;
	typedef typename utility::vector1< T >::const_iterator iter;

	assert( n_wanted <= vec.size() );

	utility::vector1< T > values;
	values.reserve( n_wanted );
	Size n_seen(0);
	for ( iter it = vec.begin(), end = vec.end(); it != end; ++it ) {
		++n_seen;
		if ( values.size() < n_wanted ) {
			values.push_back( *it );
		} else {
			Real const accept_prob(
				reservoir_sample_accept_prob( n_wanted, n_seen )
			);

			if ( rg.uniform() <= accept_prob ) {
				Size const idx = rg.random_range( 1, values.size() );
				values[ idx ] = *it;
			}
		}
	}
	return values;
} // reservoir_sample

} // namespace random
} // namespace numeric

#endif // INCLUDED_numeric_random_reservoir_sampling_HH
