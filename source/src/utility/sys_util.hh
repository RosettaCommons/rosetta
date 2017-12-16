// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/sys_util.hh
/// @brief  All system functions in utility that have no other home
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_utility_sys_util_hh
#define INCLUDED_utility_sys_util_hh


// Package headers
//#include <utility/exit.hh> // For historic reasons (exit used to be here)

// C++ headers
#include <iosfwd>
#include <string>

namespace utility {


/// @brief Sleep for a specified number of seconds
void
sys_sleep( double const seconds );

void
rand_sleep();

/// @brief Generate timestamp string
std::string
timestamp();

/// @brief Generate timestamp string, short format
std::string
timestamp_short();

} // namespace utility


#endif // INCLUDED_utility_sys_util_HH
