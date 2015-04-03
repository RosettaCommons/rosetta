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

#ifdef WIN_PYROSETTA
	#include <basetsd.h>
	/// @brief Fixed size types
	typedef  INT32  int32_t; // 32-bit unsigned integer
	typedef  UINT32 uint32_t; // 32-bit unsigned integer
	typedef  INT64  int64_t; // 64-bit signed integer
	typedef  UINT64  uint64_t; // 64-bit unsigned integer
#endif

/// @brief Fixed size types:
// int64_t  64-bit signed integer
// uint64_t  64-bit unsigned integer


/// @brief Scalable size types
// intptr_t  Pointer-sized signed integer
// uintptr_t  Pointer-sized unsigned integer
//typedef  long int  ssize_t; // Signed size

// #ifdef WIN_PYROSETTA
// 	#ifndef HAVE_SSIZE_T
// 		#define HAVE_SSIZE_T
// 		typedef  int  ssize_t; // Signed size
// 	#endif
// #endif

// #ifndef WIN_PYROSETTA
// typedef  int  ssize_t;
// #endif


// Windows native version
// #include <cstddef>
// #if defined(_MSC_VER)
// 	#include <basetsd.h>
// 	typedef SSIZE_T ssize_t;
// #endif

namespace platform {

#ifdef WIN_PYROSETTA
	// conservative, cross-compile version, we can't really use size_t here because it definition is vary between host and target compilers
	typedef unsigned int Size;
	typedef int          SSize;
	typedef unsigned int uint;


    // alsmost size_t, we can't use this version because GCC.64Bit and MSVC.64Bit defines size_t differently
	// typedef std::size_t  Size;
	// typedef ssize_t      SSize;
	// typedef std::size_t  uint;
	//unsigned __int64

#else
	typedef unsigned int  Size;
	typedef int           SSize;
	typedef unsigned int  uint;
#endif

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
