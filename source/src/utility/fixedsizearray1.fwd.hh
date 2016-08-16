// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  utility/fixedsizearray.fwd.hh
/// @brief unresizable vector whose size is known at compile time,
/// which may be allocated on the stack, and which indexes from 1.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_fixedsizearray1_fwd_hh
#define INCLUDED_utility_fixedsizearray1_fwd_hh

#include <platform/types.hh>

namespace utility {

template< typename T, platform::Size S >
class fixedsizearray1iterator;

template< typename T, platform::Size S >
class fixedsizearray1const_iterator;

template < typename T, platform::Size S >
class fixedsizearray1;


} // namespace utility

#endif
