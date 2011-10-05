// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/prob_util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef NUMERIC_PROB_UTIL_hh_
#define NUMERIC_PROB_UTIL_hh_

// C/C++ headers
#include <algorithm>
#include <iterator>

/// A collection of functions for working with probabilities
namespace numeric {

/// @brief Returns the sum of all elements on the range [first, last)
template <class InputIterator>
double sum(InputIterator first, InputIterator last) {
  double sum = 0;
  for (; first != last; ++first)
    sum += *first;
  return sum;
}

/// @brief Normalizes elements on the range [first, last)
template <class InputIterator>
void normalize(InputIterator first, InputIterator last) {
  const double div = sum(first, last);
  for (; first != last; ++first) {
    *first /= div;
  }
}

}  // namespace numeric

#endif  // NUMERIC_PROB_UTIL_hh_
