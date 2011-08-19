// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/util.hh
/// @brief small bundle of utilities for dealing with numbers.
/// @author James Thompson

#ifndef INCLUDED_numeric_util_hh
#define INCLUDED_numeric_util_hh


#include <numeric/types.hh>
#include <utility/vector1.hh>

#include <limits>
#include <cmath>
#include <algorithm>

namespace numeric {

  /// @brief Computes log(x) in the given base
  inline double log(double x, double base) {
		return std::log10(x) / std::log10(base);
  }

	/// @brief portable check to see if a value is NaN.
	template < typename T >
	inline bool isnan( T value ) {
		return value != value;
	}

	template < typename T >
	inline bool isinf( T value ) {
		return std::numeric_limits< T >::has_infinity &&
			value == std::numeric_limits< T >::infinity();
	}

	/// @brief are two Real values are equal up to some epsilon
	///
	/// implemented only for Reals, to prevent unsigned hassle
	/// (Barak 30/6/2009)
	inline bool equal_by_epsilon(
		numeric::Real value1,
		numeric::Real value2,
		numeric::Real epsilon
	) {
		if(epsilon < 0)
			epsilon = -epsilon;
		return (
			value1 <= value2 + epsilon &&
			value1 >= value2 - epsilon
		);
	}


	/// @brief Returns the median from a vector1 of Real values.
	numeric::Real median( utility::vector1< numeric::Real > const & values );

	numeric::Real mean( utility::vector1< numeric::Real > const & values );

	/// @brief Calculates the acceptance probability of a given score-change at
	/// the given temperature, generally used in simulated annealing algorithms.
	/// Returns a value in the range (0-1).
	inline
	Real boltzmann_accept_probability(
		Real const score_before,
		Real const score_after,
		Real const temperature
	) {
		Real const boltz_factor( ( score_before - score_after ) / temperature );
		Real const probability (
			std::exp( std::min( 40.0, std::max( -40.0, boltz_factor ) ) )
		);
		return probability;
	}
} // numeric

#endif
