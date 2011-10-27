// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   platform/windows/32/gcc/platform/types.hh
/// @brief  Platform-specific types for 32-bit Windows with GCC compiling PyRosetta source generated on 64-bit Linux
/// @author Sergey Lyskov


#ifndef INCLUDED_platform_pyrosetta_windows_32_gcc_platform_types_hh
#define INCLUDED_platform_pyrosetta_windows_32_gcc_platform_types_hh


// C++ headers
//#include <stdint.h> // int64_t, uint64_t, intptr_t, uintptr_t
#include <basetsd.h>

/// @brief Fixed size types
typedef  INT32  int32_t; // 32-bit unsigned integer
typedef  UINT32 uint32_t; // 32-bit unsigned integer
typedef  INT64  int64_t; // 64-bit signed integer
typedef  UINT64  uint64_t; // 64-bit unsigned integer


/// @brief Fixed size types:
// int64_t  64-bit signed integer
// uint64_t  64-bit unsigned integer


/// @brief Scalable size types
// intptr_t  Pointer-sized signed integer
// uintptr_t  Pointer-sized unsigned integer
//typedef  long int  ssize_t; // Signed size

#ifndef HAVE_SSIZE_T
    #define HAVE_SSIZE_T
	typedef  long int  ssize_t; // Signed size
#endif

namespace platform {

typedef size_t       Size;
typedef long int     SSize;
typedef size_t       uint;

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
