// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/minmax.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_UTILITY_MINMAX_HH
#define INCLUDED_UTILITY_MINMAX_HH

// C/C++ headers
#include <vector>

// Package headers
#include <utility/vector1.hh>

namespace utility {

// -- argmin -- //

/// @brief Returns the argument whose value is minimal according to operator<.
/// @details Adheres to STL numbering (0-indexed).
template <typename T>
int argmin(const std::vector<T>& iterable) {
	if ( iterable.size() == 0 ) {
		return -1;
	}

	int min_idx = 0;
	for ( std::size_t curr_idx = 1; curr_idx < iterable.size(); ++curr_idx ) {
		if ( iterable[curr_idx] < iterable[min_idx] ) {
			min_idx = curr_idx;
		}
	}

	return min_idx;
}

/// @brief Returns the argument whose value is minimal according to operator<.
/// @details Adheres to Rosetta numbering (1-indexed).
template <typename T>
int argmin(const utility::vector1<T>& iterable) {
	if ( iterable.size() == 0 ) {
		return 0;
	}

	int min_idx = 1;
	for ( std::size_t curr_idx = 2; curr_idx <= iterable.size(); ++curr_idx ) {
		if ( iterable[curr_idx] < iterable[min_idx] ) {
			min_idx = curr_idx;
		}
	}

	return min_idx;
}

// -- argmax -- //

/// @brief Returns the argument whose value is maximal according to operator>.
/// @details Adheres to STL numbering (0-indexed).
template <typename T>
int argmax(const std::vector<T>& iterable) {
	if ( iterable.size() == 0 ) {
		return -1;
	}

	int max_idx = 0;
	for ( std::size_t curr_idx = 1; curr_idx < iterable.size(); ++curr_idx ) {
		if ( iterable[curr_idx] > iterable[max_idx] ) {
			max_idx = curr_idx;
		}
	}

	return max_idx;
}

/// @brief Returns the argument whose value is maximal according to operator>.
/// @details Adheres to Rosetta numbering (1-indexed).
template <typename T>
int argmax(const utility::vector1<T>& iterable) {
	if ( iterable.size() == 0 ) {
		return 0;
	}

	int max_idx = 1;
	for ( std::size_t curr_idx = 2; curr_idx <= iterable.size(); ++curr_idx ) {
		if ( iterable[curr_idx] > iterable[max_idx] ) {
			max_idx = curr_idx;
		}
	}

	return max_idx;
}

}  // namespace utility

#endif  // INCLUDED_UTILITY_MINMAX_HH
