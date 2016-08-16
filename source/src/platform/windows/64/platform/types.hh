// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   platform/windows/64/platform/types.hh
/// @brief  Platform-specific types for 64-bit Windows
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_platform_windows_64_platform_types_hh
#define INCLUDED_platform_windows_64_platform_types_hh


/// @brief Fixed size types
typedef  signed int int32_t; // 32-bit signed integer
typedef  unsigned int uint32_t; // 32-bit unsigned integer
typedef  long long int  int64_t; // 64-bit signed integer
typedef  unsigned long long int  uint64_t; // 64-bit unsigned integer


/// @brief Scalable size types
typedef  long long int  intptr_t; // Pointer-sized signed integer
typedef  unsigned long long int  uintptr_t; // Pointer-sized unsigned integer
typedef  long long int  ssize_t; // Signed size


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
bool const CASE_SENSITIVE( false );

/// @brief Volume specifier used?
bool const VOLUME_USED( true );

/// @brief File path separator
char const PATH_SEPARATOR( '\\' );


} // namespace file


} // namespace platform


#endif // INCLUDED_platform_types_HH
