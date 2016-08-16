// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   platform/windows/32/gcc/platform/types.hh
/// @brief  Platform-specific types for 32-bit Windows with GCC
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_platform_windows_32_gcc_platform_types_hh
#define INCLUDED_platform_windows_32_gcc_platform_types_hh


// C++ headers
#include <stdint.h> // int64_t, uint64_t, intptr_t, uintptr_t


/// @brief Fixed size types:
// int64_t  64-bit signed integer
// uint64_t  64-bit unsigned integer


/// @brief Scalable size types
// intptr_t  Pointer-sized signed integer
// uintptr_t  Pointer-sized unsigned integer
typedef  long int  ssize_t; // Signed size


namespace platform {


namespace file {


/// @brief Are file/path names case sensitive?
bool const CASE_SENSITIVE( false );

/// @brief Volume specifier used?
bool const VOLUME_USED( true );

/// @brief File path separator
char const PATH_SEPARATOR( '\\' );


} // namespace file


} // namespace platform


#endif // INCLUDED_platform_types_HH
