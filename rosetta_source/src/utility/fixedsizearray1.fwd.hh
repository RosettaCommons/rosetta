// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
