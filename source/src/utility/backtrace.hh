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

// Provide workaround for modern compiler feature, if missing
#ifdef __has_include
#define MY__has_include( x ) __has_include( x )
#else
#define MY__has_include( x ) 1
#endif

// Shared with utility/exit.hh
#ifndef NORETURN
#ifdef __GNUC__
#  define NORETURN __attribute__ ((noreturn))
#elif __clang__
#  define NORETURN __attribute__ ((noreturn))
#else
#  define NORETURN
#endif
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

#include <execinfo.h>
#include <cxxabi.h>
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <iostream>

#include <utility/CSI_Sequence.hh>


//////////////////////////////////////////////////////////////////////
//
// See note above if this is causing your build to not compile.
//
// from user 'john' in stackoverflow, but adapted to work on a Mac.
//  and then also re-adapted to work on linux.
//  -- rhiju, 2014
//////////////////////////////////////////////////////////////////////
inline
std::string
demangle( std::string trace ) {

	std::string::size_type begin, end;

	// find the beginning and the end of the useful part of the trace

	// On a Mac, part to be demangled starts with underscore, and then is followed by "+" increment,
	//  separated by spaces from surrounding.
	begin = trace.find(" _") + 1;
	end   = trace.find(" +",begin);

	//How it looks for Linux, with parentheses around part to be demangled.
	// /home/rhiju/src/rosetta/main/source/cmake/build_release/libutility.so(_Z15print_backtracev+0x23) [0x7fb75e08c1a3]
	if ( begin == std::string::npos || end == std::string::npos ) {
		begin = trace.find("(_") + 1;
		end   = trace.find("+",begin);
	}

	// if begina and end were found, we'll go ahead and demangle
	if ( begin != std::string::npos && end != std::string::npos ) {
		std::string mangled_trace = trace.substr(begin, end - begin);
		size_t maxName = 1024;
		int demangleStatus;

		// If output_buffer to __cxa_demangle() is null, memory will be allocated and the pointer returned.
		char* demangledName; // = (char*) malloc(maxName);
		if ( (demangledName = abi::__cxa_demangle(mangled_trace.c_str(), 0, &maxName,
				&demangleStatus)) && demangleStatus == 0 ) {
			trace = trace.substr(0,begin) + demangledName + trace.substr(end ); // the demangled name is now in our trace string
		}
		free(demangledName); // Will handle null pointers gracefully.
	}
	return trace;
}

////////////////////////////////////////////////////////////////////
//
// See note above if this is causing your build to not compile.
//
// stolen directly from
//
// https://developer.apple.com/library/mac/documentation/Darwin/Reference/ManPages/man3/backtrace.3.html
//
// see also:
//
// http://man7.org/linux/man-pages/man3/backtrace.3.html#NOTES
//
//   -- rhiju, 2014.
////////////////////////////////////////////////////////////////////

inline
bool
print_backtrace( char const * condition ) {

	// instead of printing a backtrace and letting the assert() function exectue,
	// (and halt the program) throw an exception.  Look, don't rely on this
	// functionality in anything besides your unit tests.
	maybe_throw_on_next_assertion_failure( condition );

	size_t const callstack_size = 128;
	void* callstack[callstack_size];
	int i, frames = backtrace(callstack, callstack_size);
	char** strs = backtrace_symbols(callstack, frames);
	std::cerr << utility::CSI_Magenta(); // set color of cerr to magenta
	for ( i = 0; i < frames; ++i ) {
		std::cerr << demangle( strs[i] ).c_str() << std::endl;
	}
	std::cerr << utility::CSI_Reset(); // reset color of cerr
	free(strs);
	return false; // allows use in debug_assert
}

/// @brief A version of print_backtrace specifically for debug_assert, which
/// tells the Clang Static Analyzer that we shouldn't continue on if we hit this point.
inline bool print_backtrace_NR( char const * condition ) NORETURN;

inline
bool
print_backtrace_NR( char const * condition ) {
	print_backtrace( condition );
	assert(false);
	abort(); // To make the compiler happy on release-mode builds
}

#define debug_assert(condition) {assert( ( condition ) || print_backtrace_NR( #condition ) ); }

#else
// _WIN32, etc.
#include <assert.h>

inline
void
print_backtrace( char const * /*unnamed*/ ){
	// no op
	// if someone cares, should be possible to code up a backtrace for Windows!
	// function signature needs to match windows and *nix builds.
}

#define debug_assert(condition) {assert( condition || maybe_throw_on_next_assertion_failure( #condition ) ); }

#endif

#endif // INCLUDED_utility_backtrace_HH
