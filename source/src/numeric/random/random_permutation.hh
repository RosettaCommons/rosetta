// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/random_permutation.hh
/// @brief  Reorder the elements in a vector using a RNG object
/// @author Andrew Leaver-Fay
///
/// @remarks
///  @li -


#ifndef INCLUDED_numeric_random_random_permutation_hh
#define INCLUDED_numeric_random_random_permutation_hh


// Unit headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <algorithm>
#include <cmath>
#include <iostream>

namespace numeric {
namespace random {

template< class T >
void
random_permutation(
	utility::vector1< T > & vect,
	RandomGenerator & rg
)
{
	Size const vsize = vect.size();
	for ( Size ii = vsize; ii >= 2; --ii ) {
		Size swap_partner = static_cast< Size > (std::floor( ii*rg.uniform() + 1 ));

		T temp = vect[ swap_partner ];
		vect[ swap_partner ] = vect[ ii ];
		vect[ ii ] = temp;
	}
}

template< class T >
void
random_permutation(
	utility::vector0< T > & vect,
	RandomGenerator & rg
)
{
	Size const vsize = vect.size();
	for ( Size ii = vsize - 1; ii >= 1; --ii ) {
		Size swap_partner = static_cast< Size > ( std::floor( ii*rg.uniform() ));

		T temp = vect[ swap_partner ];
		vect[ swap_partner ] = vect[ ii ];
		vect[ ii ] = temp;
	}
}


template< class T >
void
random_permutation(
	std::vector< T > & vect,
	RandomGenerator & rg
)
{
	Size const vsize = vect.size();
	for ( Size ii = vsize - 1; ii >= 1; --ii ) {
		Size swap_partner = static_cast< Size > ( std::floor( ii*rg.uniform() ));

		T temp = vect[ swap_partner ];
		vect[ swap_partner ] = vect[ ii ];
		vect[ ii ] = temp;
	}
}


/// @brief randomly shuffle elements of a sequence
template< typename RandomAccessIterator >
void
random_permutation(
	RandomAccessIterator first,
	RandomAccessIterator last,
	RandomGenerator & rg
)
{
	if ( first != last ) {
		for ( RandomAccessIterator i = first + 1; i != last; ++i ) {
			std::iter_swap( i, first + static_cast< Size >( ( i - first + 1 ) * rg.uniform() ) );
		}
	}
}


} // namespace random
} // namespace numeric


#endif // INCLUDED_numeric_random_random_permutation_HH
