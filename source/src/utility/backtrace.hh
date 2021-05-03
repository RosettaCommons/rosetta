// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/backtrace.hh
/// @brief  Programmatic backtrace whenever you want it.
/// @author Rhiju Das

#ifndef INCLUDED_utility_backtrace_hh
#define INCLUDED_utility_backtrace_hh

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Note to devs -- if your build is failing to compile because
// of unrecognized functions backtrace() or backtrace_symbols(),
// then just expand the
//
//   #ifdef _WIN32
//
// to include some tag that goes with your build.
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

#include <cassert>
#include <string>

// Provide workaround for modern compiler feature, if missing
#ifdef __has_include
#define MY__has_include( x ) __has_include( x )
#else
#define MY__has_include( x ) 1
#endif

/// @brief Function for unit testing only -- if an assertion failure is hit, throw an exception
/// instead of exiting.  Don't let me catch you calling this function from anywhere besides a
/// unit test.  Punishment will be swift.
void set_throw_on_next_assertion_failure();

/// @brief Throw an exception if set_throw_on_next_assertion_failure was called since the
/// last time this function was called.
bool maybe_throw_on_next_assertion_failure( char const * condition );


// C++ headers
#if defined(__GNUC__)  &&  !defined(WIN32)  &&  !defined(__CYGWIN__) && MY__has_include( <cxxabi.h> ) && !defined(ANDROID)

std::string
demangle( std::string trace );

std::string
backtrace_string(int skip=0);

bool
print_backtrace( char const * /*unused*/ );

#else
// _WIN32, etc.
#include <assert.h>

inline
std::string
backtrace_string(int skip=0) {
	return "";
}

inline
bool
print_backtrace( char const * /*unused*/ ) {
	// no op
	// if someone cares, should be possible to code up a backtrace for Windows!
	// function signature needs to match windows and *nix builds.
	return false; // allows use in debug_assert
}

#endif // Platform testing

#ifdef __clang_analyzer__
/// @details NORETURN_ATTR to tell the Clang Static Analyzer that we don't continue on if we hit this point.
NORETURN_ATTR
#endif
bool handle_assert_failure( char const * condition, std::string const & file, int const line );

// When used, this macro must be followed by a semi-colon to be beautified properly.
#define debug_assert(condition) assert( ( condition ) || handle_assert_failure( #condition, __FILE__, __LINE__ ) )


#endif // INCLUDED_utility_backtrace_HH
