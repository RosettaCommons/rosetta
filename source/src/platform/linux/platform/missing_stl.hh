// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   platform/linux/platform/missing_stl.hh
/// @brief  Missing STL functions
/// @author David E. Kim (dekim@u.washington.edu)

/*
The Android NDK (version android-ndk-r14b) gnustl is missing some functions so include them
here until Google adds them. The previous solution was to use Crystax which has better c++11
support but it crashes with Andoid 7 (nougat). The Android NDK LLVM libc++ also had issues.

Info about Android NDK C++ library support can be found at:

https://developer.android.com/ndk/guides/cpp-support.html

*/

#ifndef INCLUDED_platform_linux_platform_missing_stl_hh
#define INCLUDED_platform_linux_platform_missing_stl_hh

#ifdef ANDROID





#include <tgmath.h>
#include <stdlib.h>
#include <utility/string_util.hh>

namespace std
{

template <class T>
inline platform::Real
round(const T & r)
{
	return ::round(r);
}

inline platform::Size
stoi(const string & s)
{
	return ::atoi(s.c_str());
}

inline platform::Real
stod(const string & s)
{
	return ::atof(s.c_str());
}

inline platform::Real
stof(const string & s)
{
	return ::atof(s.c_str());
}

template <class T>
inline string
to_string (const T & t)
{
	return utility::to_string(t);
}

}

#endif

#endif // INCLUDED_platform_missing_stl_HH
