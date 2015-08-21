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
#include <numeric/numeric.functions.hh>
#include <utility/vector1.hh>

#include <limits>
#include <cmath>
#include <algorithm>

namespace numeric {

/// @brief Clamps <value> to the closed interval [lower_bound, upper_bound].
/// Templated type must implement operator<.
template<typename Number>
Number clamp(Number value, Number lower_bound, Number upper_bound) {
	if ( value < lower_bound ) {
		return lower_bound;
	} else if ( upper_bound < value ) {
		return upper_bound;
	} else {
		return value;
	}
}

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
	if ( epsilon < 0 ) {
		epsilon = -epsilon;
	}
	return (
		value1 <= value2 + epsilon &&
		value1 >= value2 - epsilon
	);
}


/// @brief Returns the median from a vector1 of Real values.
numeric::Real median( utility::vector1< numeric::Real > const & values );

numeric::Real mean( utility::vector1< numeric::Real > const & values );

template < typename T >
T max( utility::vector1< T > const & values ) {
	assert(values.size());
	T retval = values[1];
	for ( numeric::Size ii(2); ii <= values.size(); ++ii ) {
		retval = max(retval, values[ii]);
	}
	return retval;
}

template < typename T >
T min( utility::vector1< T > const & values ) {
	assert(values.size());
	T retval = values[1];
	for ( numeric::Size ii(2); ii <= values.size(); ++ii ) {
		retval = min(retval, values[ii]);
	}
	return retval;
}


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

/// @brief recursive binary search that finds the value closest to key.  Call find_nearest_value(input_list,value) instead.  It's the driver
/// function for this function.  This fails miserably (and silently!) on a non-sorted vector, so don't do that!.
template < typename T >
T find_nearest_value(typename utility::vector1<T> const & input_list, T key, platform::Size min_index, platform::Size max_index)
{

	if ( (max_index - min_index) == 1 ) {
		//down to the last two items and we haven't found the key yet, one of these must be the closest.
		T max_distance = std::abs(input_list[max_index] - key);
		T min_distance = std::abs(input_list[min_index] - key);
		if ( max_distance < min_distance ) {
			return input_list[max_index];
		} else {
			return input_list[min_index];
		}
	}

	platform::Size mid_index = (min_index + max_index) /2;

	if ( input_list[mid_index] > key ) {
		return find_nearest_value<T>(input_list,key,min_index,mid_index-1);
	} else if ( input_list[mid_index] < key ) {
		return find_nearest_value<T>(input_list,key,mid_index+1, max_index);
	} else {
		return input_list[mid_index];
	}
}

/// @brief given a vector and an input value, return the value in the vector that is closest to the input
/// This is a wrapper for find_nearest_value(input_list,key,min,max) and insures that you're sorted.
template < typename T >
T find_nearest_value(typename utility::vector1<T> const & input_list, T key)
{
	utility::vector1<T> temp_list = input_list;
	std::sort(temp_list.begin(),temp_list.end());

	platform::Size initial_minimum = 1;
	platform::Size initial_maximum = temp_list.size();

	return find_nearest_value<T>(temp_list,key,initial_minimum,initial_maximum);

}

} // numeric

#endif
