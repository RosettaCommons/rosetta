// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/vector1.functions.hh
/// @brief  some useful functions for vector1 that I'm tired of rewriting
/// @author Andrew Leaver-Fay ( aleaverfay@gmail.com )


#ifndef INCLUDED_utility_vector1_functions_HH
#define INCLUDED_utility_vector1_functions_HH


// Unit headers
#include <utility/vector1.hh>

#include <utility/numbers.hh>

namespace utility {

/// @brief Find the largest value in a vector
///
/// class T must provide an operator < ()  and operator = ().
/// Error if input.size() == 0
template < class T >
T
max( vector1< T > const & input )
{
	debug_assert( input.size() > 0 );

	T largest_so_far = input[ 1 ];
	for ( typename vector1< T >::Size ii = 2; ii <= input.size(); ++ii ) {
		if ( largest_so_far < input[ ii ] ) largest_so_far = input[ ii ];
	}
	return largest_so_far;
}

/// @brief Find the smallest value in a vector
///
/// class T must provide an operator < ()  and operator = ().
/// Error if input.size() == 0
template < class T >
T
min( vector1< T > const & input )
{
	debug_assert( input.size() > 0 );

	T smallest_so_far = input[ 1 ];
	for ( typename vector1< T >::Size ii = 2; ii <= input.size(); ++ii ) {
		if ( input[ ii ] < smallest_so_far ) smallest_so_far = input[ ii ];
	}
	return smallest_so_far;
}

/// @brief find the index of the largest value in a vector
///
/// class T must provide an operator < ()  and operator = ().
/// Error if input.size() == 0
template < class T >
typename vector1< T >::Size
arg_max( vector1< T > const & input )
{
	debug_assert( input.size() > 0 );

	T largest_so_far = input[ 1 ];
	typename vector1< T >::Size index_of_largest = 1;
	for ( typename vector1< T >::Size ii = 2; ii <= input.size(); ++ii ) {
		if ( largest_so_far < input[ ii ] ) {
			largest_so_far = input[ ii ];
			index_of_largest = ii;
		}
	}
	return index_of_largest;
}

/// @brief find the index of the smallest value in a vector
///
/// class T must provide an operator < ()  and operator = ().
/// Error if input.size() == 0
template < class T >
typename vector1< T >::Size
arg_min( vector1< T > const & input )
{
	debug_assert( input.size() > 0 );

	T smallest_so_far = input[ 1 ];
	typename vector1< T >::Size index_of_smallest = 1;
	for ( typename vector1< T >::Size ii = 2; ii <= input.size(); ++ii ) {
		if ( input[ ii ] < smallest_so_far ) {
			smallest_so_far = input[ ii ];
			index_of_smallest = ii;
		}
	}
	return index_of_smallest;
}

template < class T >
void
insert_middle(
	vector1< T > & vect,
	typename vector1< T >::Size const position,
	T const new_value,
	bool expand
)
{
	if ( position == vect.size() + 1 ) {
		if ( expand ) vect.push_back( new_value );
	} else if ( position == vect.size() ) {
		if ( expand ) {
			vect.push_back( vect[ vect.size() ] );
			vect[ position ] = new_value;
		}
	} else {
		T next_val = vect[ position ];
		vect[ position ] = new_value;
		for ( typename vector1< T >::Size ii = position+1; ii <= vect.size(); ++ii ) {
			std::swap( vect[ ii ], next_val );
		}
		if ( expand ) {
			vect.push_back( next_val );
		}
	}
}

/// @brief Finds indices of the n largest items in input vector, where
/// n is size of the arg_list vector.  The indices are reported in decreasing sorted
/// order by the value of the corresponding position in the input vector.
/// If m is the size of the input vector, then this function is O( m n**2)
///
/// @details If arg_list has more entries than input, then it is resized
/// at the end of this function to have as many entries as the input.
template < class T >
void
arg_greatest_several(
	vector1< T > const & input,
	vector1< typename vector1< T >::Size > & arg_list )
{
	typename vector1< T >::Size const n = arg_list.size();
	arg_list.resize( 0 );
	for ( typename vector1< T >::Size ii = 1; ii <= input.size(); ++ii ) {
		if ( arg_list.size() == 0 && n != 0 ) {
			arg_list.push_back( ii );
		} else {
			for ( typename vector1< T >::Size jj = arg_list.size(); jj >= 1; --jj ) {
				if ( input[ ii ] <= input[ arg_list[ jj ] ] ) {
					insert_middle( arg_list, jj + 1, ii, arg_list.size() < n );
					break;
				} else if ( jj == 1 ) {
					insert_middle( arg_list, 1, ii, arg_list.size() < n );
				}
			}
		}
	}
}


/// @brief Finds indices of the n smallest items in input vector, where
/// n is size of the arg_list vector.  The indices are reported in increasing sorted
/// order by the value of the corresponding position in the input vector.
/// If m is the size of the input vector, then this function is O( m n**2)
///
/// @details If arg_list has more entries than input, then it is resized
/// at the end of this function to have as many entries as the input.
template < class T >
void
arg_least_several(
	vector1< T > const & input,
	vector1< typename vector1< T >::Size > & arg_list )
{
	typename vector1< T >::Size const n = arg_list.size();
	arg_list.resize( 0 );
	for ( typename vector1< T >::Size ii = 1; ii <= input.size(); ++ii ) {
		if ( arg_list.size() == 0 ) {
			arg_list.push_back( ii );
		} else {
			for ( typename vector1< T >::Size jj = arg_list.size(); jj >= 1; --jj ) {
				if ( input[ ii ] >= input[ arg_list[ jj ] ] ) {
					insert_middle( arg_list, jj + 1, ii, arg_list.size() < n );
					break;
				} else if ( jj == 1 ) {
					insert_middle( arg_list, 1, ii, arg_list.size() < n );
				}
			}
		}
	}
}

//template< class T >
//void
//erase(vector1< T > & input, core::Size index)
//{
//
// input::iterator iter_this = std::find( input.begin(), input.end(), input[ index ] );
// assert( iter_this != input.end() );
// input.erase( iter_this );
//
//}


/// @brief Find the index into the input range array such that
/// range_array[ index ] <= query < range_array[ index + 1 ] or if
/// index+1 would be off the end of the array, then simply the last
/// index in the array (for which range_array[ index ] <= query holds).
/// Each entry in range_array represents the lower boundry (inclusive)
/// for a range, so that range i is defined by
/// ( range_array[ i ], range_array[i+1] - 1 ) except for the last
/// range, which is from range_array[ range_array.size() ] to infinity.
/// T must be a discrete type.
/// range_array must be sorted. range_array[ 1 ] must be as small as
/// the smallest possible query.
template < typename T >
platform::Size
binary_search_ranges(
	utility::vector1< T > const & range_array,
	T query
)
{
	platform::Size lower( 1 ), upper( range_array.size() );
	platform::Size guess( 0 );

	while ( true ) {

		guess = ( upper + lower ) / 2;
		if ( guess == upper ) return upper;

		T guess_val = range_array[ guess ];
		T next_val = range_array[ guess+1 ];
		if ( guess_val <= query && query < next_val ) {
			// found it!
			break;
		}
		if ( guess_val <= query ) {
			lower = guess+1;
		} else {
			upper = guess-1;
		}
	}

	return guess;

}

/// @brief Are any of the values contained in the vector non-finite (NaN, inf)
template< platform::SSize L, typename T >
bool
is_finite( utility::vectorL< L, T > const & vect )
{
	for ( T const & v: vect ) {
		if ( ! utility::isfinite(v) ) {
			return false;
		}
	}
	return true;
}


} // namespace utility

#endif
