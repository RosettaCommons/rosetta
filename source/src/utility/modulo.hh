// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/modulo.hh
/// @brief  calculates modulo for integer, since C++ '%' operator has weird behavior for negative numbers.
/// @author Rhiju Das


#ifndef INCLUDED_utility_modulo_hh
#define INCLUDED_utility_modulo_hh

// Unit headers
#include <platform/types.hh>

// C++ headers
#include <stdlib.h>

namespace utility {

/// @brief modulo of an input integer a (can be negative) with respect to unsigned integer b

/// @details
///   Most folks use C's "a % b" operator but it gives funny behavior for negative a.
///   This version came out of stack overflow.
inline
platform::Size
modulo( int const & a, int const & b)
{
	return ( a >= 0 ) ? a % b : ( b - abs ( a % b ) ) % b;
}


} // namespace utility

#endif // INCLUDED_utility_modulo_HH
