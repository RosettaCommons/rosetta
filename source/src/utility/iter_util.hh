// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/iter_util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_UTILITY_ITER_UTIL_hh
#define INCLUDED_UTILITY_ITER_UTIL_hh

// C/C++ headers
#include <algorithm>
#include <iterator>

namespace utility {

/// @brief Returns an iterator on the sorted range [first, last) nearest
/// to value. If value is equidistant between adjacent elements, the
/// lesser is returned.
template <typename BidirectionalIterator, typename T>
BidirectionalIterator find_closest(BidirectionalIterator first,
	BidirectionalIterator last,
	const T& value) {
	BidirectionalIterator before = std::lower_bound(first, last, value);

	if ( before == first ) return first;
	if ( before == last )  return --last; // iterator must be bidirectional

	BidirectionalIterator after = before;
	--before;

	return (*after - value) < (value - *before) ? after : before;
}

}  // namespace utility

#endif  // INCLUDED_UTILITY_ITER_UTIL_hh
