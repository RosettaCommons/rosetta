// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/backtrace.hh
/// @brief  Programmatic backtrace whenever you want it.
/// @author Rhiju Das
///


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

// C++ headers
#if defined(__GNUC__)  &&  !defined(WIN32)  &&  !defined(__CYGWIN__)

#include <execinfo.h>
#include <cxxabi.h>
#include <string>
#include <cstdio>
#include <stdlib.h>
#include <iostream>

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
	//	/home/rhiju/src/rosetta/main/source/cmake/build_release/libutility.so(_Z15print_backtracev+0x23) [0x7fb75e08c1a3]
 	if (begin == std::string::npos || end == std::string::npos) {
		begin = trace.find("(_") + 1;
		end   = trace.find("+",begin);
	}

 	// if begina and end were found, we'll go ahead and demangle
 	if (begin != std::string::npos && end != std::string::npos) {
		std::string mangled_trace = trace.substr(begin, end - begin);
 		size_t maxName = 1024;
 		int demangleStatus;

 		char* demangledName = (char*) malloc(maxName);
 		if ((demangledName = abi::__cxa_demangle(mangled_trace.c_str(), demangledName, &maxName,
 																						 &demangleStatus)) && demangleStatus == 0) {
 			trace = trace.substr(0,begin) + demangledName + trace.substr(end ); // the demangled name is now in our trace string
 		}
		free(demangledName);
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
print_backtrace() {

  void* callstack[128];
	int i, frames = backtrace(callstack, 128);
	char** strs = backtrace_symbols(callstack, frames);
	for (i = 0; i < frames; ++i) {
		std::cerr << demangle( strs[i] ).c_str() << std::endl;
	}
	free(strs);
	return false; // allows use in debug_assert
}

#define debug_assert(condition) { assert( ( condition ) || print_backtrace() ); }

#else
// _WIN32, etc.

inline
void
print_backtrace(){
	// no op
	// if someone cares, should be possible to code up a backtrace for Windows!
}

#define debug_assert(condition) { assert( condition ); }

#endif

#endif // INCLUDED_utility_backtrace_HH
