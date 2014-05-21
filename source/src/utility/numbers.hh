// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/numbers.hh
/// @brief Brief utility classes for numeric usage
/// @details For more complex numeric utilties, see the numeric namespace
/// @author Rocco Moretti

#ifndef INCLUDED_utility_numbers_hh
#define INCLUDED_utility_numbers_hh

#include <platform/types.hh>

#include <limits>
#include <cmath>

namespace utility {

	/// @brief Get a numeric value for Size that represents an "undefined" value
	inline
	platform::Size get_undefined_size() {
		return std::numeric_limits< platform::Size >::max(); // Choice of value same as the BCL (Meiler Lab)
	}

	/// @brief Check if a Size is undefined (i.e has the same value as utility::get_undefined_size() )
	inline
	bool is_undefined( platform::Size const & val) {
		return val == get_undefined_size();
	}

	/// @brief Get a numeric value for Real that represents an "undefined" value
	inline
	platform::Real get_undefined_real() {
		return std::numeric_limits< platform::Real >::quiet_NaN(); // Choice of value same as the BCL (Meiler Lab)
	}

	/// @brief Check if a Real is undefined (i.e has the same value as utility::get_undefined_real() )
	inline
	bool is_undefined( platform::Real const & val) {
		return isnan( val ) || isinf( val );
	}


} // numeric

#endif
