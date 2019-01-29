// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/sort.functions.hh
/// @brief  Functions related to sorting
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_utility_sort_functions_hh
#define INCLUDED_utility_sort_functions_hh


#include <utility/vectorL.hh>
#include <utility/vector1.hh>

// C++ headers
#include <numeric>


namespace utility {

// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <platform::SSize L, typename T>
utility::vector1< platform::Size > argsort( utility::vectorL<L, T> const & v ) {

	// initialize original index locations
	utility::vector1< platform::Size > idx( v.size() );
	std::iota( idx.begin(), idx.end(), L );

	// sort indexes based on comparing values in v
	sort( idx.begin(), idx.end(),
		[&v]( platform::Size i1, platform::Size i2) { return v[i1] < v[i2]; } );

	return idx;
}

/// @brief Return the fractional rank of a list of numbers
template <platform::SSize L, typename T>
utility::vector1< platform::Real > fractional_rank( utility::vectorL<L, T> const & v ) {

	platform::Size total = v.size();
	utility::vector1< platform::Size > sorted_indices = argsort( v );

	// If we ignored ties: ranks[i] = sorted_indices.index( i )
	utility::vector1< platform::Real > ranks( total );

	// This part is O(n)
	for ( platform::Size i = 1; i <= total; ) {

		T const & value = v[ sorted_indices[i] ];

		// Figure out how many ties there are
		platform::Real rank_sum = i;
		platform::Size j;
		for ( j = i + 1; j <= total; j++ ) {
			if ( v[ sorted_indices[j] ] != value ) break;
			rank_sum += j;
		}
		platform::Size ties = j - i;

		// Take the average value
		rank_sum = ties == 1 ? rank_sum : rank_sum / ties;

		// Write ties into ranks
		for ( platform::Size k = i; k < j; k++ ) {
			ranks[ sorted_indices[k] ] = rank_sum;
		}

		// update i for next iteration
		i = j;
	}

	return ranks;
}


} // namespace utility


#endif // INCLUDED_utility_sort_functions_HH
