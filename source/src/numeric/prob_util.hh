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

#ifndef INCLUDED_NUMERIC_PROB_UTIL_hh
#define INCLUDED_NUMERIC_PROB_UTIL_hh

// C/C++ headers
#include <iostream>

// Utility headers
#include <utility/vector1.hh>

/// A collection of functions for working with probabilities
namespace numeric {

/// @brief Returns the sum of all elements on the range [first, last)
template <class InputIterator>
double sum(InputIterator first, InputIterator last) {
	double sum = 0;
	for ( ; first != last; ++first ) {
		sum += *first;
	}
	return sum;
}

/// @brief Normalizes elements on the range [first, last)
template <class InputIterator>
void normalize(InputIterator first, InputIterator last) {
	const double div = sum(first, last);
	for ( ; first != last; ++first ) {
		*first /= div;
	}
}

/// @brief Converts pdf to cdf
template <class RandomAccessIterator>
void cumulative(RandomAccessIterator first, RandomAccessIterator last) {
	normalize(first, last);
	for ( RandomAccessIterator i = first + 1; i != last; ++i ) {
		*i += *(i - 1);
	}
}

/// @brief Multiplies two probability vectors with one another.
/// Probability vectors are assumed to have equal lengths.
template <class ForwardIterator>
void product(ForwardIterator probs1_first, ForwardIterator probs1_last,
	ForwardIterator probs2_first, ForwardIterator probs2_last) {
	normalize(probs1_first, probs1_last);
	normalize(probs2_first, probs2_last);

	ForwardIterator i = probs1_first;
	ForwardIterator j = probs2_first;
	for ( ; i != probs1_last; ++i, ++j ) {
		*i *= *j;
	}
	normalize(probs1_first, probs1_last);
}

/// @brief Loads normalized, per-residue probabilities from filename,
/// storing the result in probs. Assumes line i holds the probability
/// of sampling residue i. There must be 1 line for each residue in the
/// pose on which this data will be used.
void read_probabilities_or_die(const std::string& filename, utility::vector1<double>* probs);

/// @brief Writes probs to the specified ostream
void print_probabilities(const utility::vector1<double>& probs, std::ostream& out);

}  // namespace numeric

#endif  // NUMERIC_PROB_UTIL_hh_
