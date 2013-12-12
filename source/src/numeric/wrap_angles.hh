// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Functions to wrap angles in different ranges.
/// author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_numeric_wrap_angles_hh
#define INCLUDED_numeric_wrap_angles_hh

#include <numeric/NumericTraits.hh>
#include <numeric/numeric.functions.hh>

namespace numeric {

/// @brief Wrap the given angle in the range [0, 2 * pi).
/// @details No conversion to radians is implied.
template<typename T>
inline T wrap_2pi(T const &angle) {
	return modulo<T>(angle, NumericTraits<T>::pi_2());
}

/// @brief Wrap the given angle in the range [-pi, pi).
/// @details No conversion to radians is implied.
template<typename T>
inline T wrap_pi(T const &angle) {
	return wrap_2pi<T>(angle + NumericTraits<T>::pi()) - NumericTraits<T>::pi();
}

/// @brief Wrap the given angle in the range [0, 360).
/// @details No conversion to degrees is implied.
template<typename T>
inline T wrap_360(T const &angle) {
	return modulo<T>(angle, 360);
}

/// @brief Wrap the given angle in the range [-180, 180).
/// @details No conversion to degrees is implied.
template<typename T>
inline T wrap_180(T const &angle) {
	return modulo<T>(angle + 180, 360) - 180;
}

} // end namespace numeric

#endif
