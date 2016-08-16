// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   platform/macos/platform/types.hh
/// @brief  Platform-specific types for MacOS 10.4 on 32-bit PPC with GCC
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


#ifndef INCLUDED_platform_macos_platform_types_hh
#define INCLUDED_platform_macos_platform_types_hh


// C++ headers
#include <stdint.h>
#include <sys/types.h>
#include <cstddef> // std::size_t

/// @brief Fixed size types
// int64_t  64-bit signed integer
// uint64_t  64-bit unsigned integer


/// @brief Scalable size types
// intptr_t  Pointer-sized signed integer
// uintptr_t  Pointer-sized unsigned integer
// ssize_t  Signed size


namespace platform {

typedef std::size_t  Size;
typedef ssize_t      SSize;
typedef std::size_t  uint;

// Floating point precision control scalar
#ifdef ROSETTA_FLOAT // Real == float
typedef  float  Real;
#else // Real == double
typedef  double  Real;
#endif

namespace file {


/// @brief Are file/path names case sensitive?
bool const CASE_SENSITIVE( true );

/// @brief Volume specifier used?
bool const VOLUME_USED( false );

/// @brief File path separator
char const PATH_SEPARATOR( '/' );


} // namespace file


} // namespace platform


#endif // INCLUDED_platform_types_HH
