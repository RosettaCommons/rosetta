// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/backtrace.cc
/// @brief  Instead of printing a backtrace inside of an assertion failure, throw
///         an exception.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <utility/backtrace.hh>
#include <utility/crash_report.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/CSI_Sequence.hh>

#include <iostream>
#include <sstream>

#if defined(__GNUC__)  &&  !defined(WIN32)  &&  !defined(__CYGWIN__) && MY__has_include( <cxxabi.h> ) && !defined(ANDROID)

#include <execinfo.h>
#include <cxxabi.h>
#include <string>
#include <cstdio>
#include <stdlib.h>

#endif

static bool throw_the_next_time_an_assertion_failure_is_hit( false );

void set_throw_on_next_assertion_failure()
{
	throw_the_next_time_an_assertion_failure_is_hit = true;
}

bool maybe_throw_on_next_assertion_failure( char const * condition )
{
	if ( throw_the_next_time_an_assertion_failure_is_hit ) {
		throw_the_next_time_an_assertion_failure_is_hit = false;
		throw CREATE_EXCEPTION(utility::excn::Exception, std::string( "assertion failure hit:" ) + condition );
	}
	return false;
}

bool
handle_assert_failure( char const * condition, std::string const & file, int const line ) {

	// Instead of printing a backtrace and halting the program, throw an exception.
	// Look, don't rely on this functionality in anything besides your unit tests.
	maybe_throw_on_next_assertion_failure( condition );

	std::ostringstream oss;
	oss << "\nERROR: Assertion `" << condition << "` failed.\n";
	oss << "ERROR:: Exit from: " << file << " line: " << line << "\n";

	std::cerr << utility::CSI_Reset() << utility::CSI_Red() << utility::CSI_Bold();
	std::cerr << oss.str();
	std::cerr << utility::CSI_Reset();
	std::cerr << std::endl;

	//print_backtrace( condition );
	utility::save_crash_report(oss.str(), file, line);


#ifdef __clang_analyzer__
	abort(); // To make the compiler happy on release-mode builds
#else
	return false;
#endif
}


#if defined(__GNUC__)  &&  !defined(WIN32)  &&  !defined(__CYGWIN__) && MY__has_include( <cxxabi.h> ) && !defined(ANDROID)

bool
print_backtrace( char const * /*unused*/ ) {
	std::cerr << utility::CSI_Magenta(); // set color of cerr to magenta
	std::cerr << "BACKTRACE:\n";
	std::cerr <<  backtrace_string();
	std::cerr << utility::CSI_Reset();
	return false; // allows use in debug_assert
}

//////////////////////////////////////////////////////////////////////
//
// See note in utility/backtrace.hh if this is causing your build to not compile.
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
		if ( (demangledName = abi::__cxa_demangle(mangled_trace.c_str(), nullptr, &maxName,
				&demangleStatus)) && demangleStatus == 0 ) {
			trace = trace.substr(0,begin) + demangledName + trace.substr(end ); // the demangled name is now in our trace string
		}
		free(demangledName); // Will handle null pointers gracefully.
	}
	return trace;
}

////////////////////////////////////////////////////////////////////
//
// See note in utility/backtrace.hh if this is causing your build to not compile.
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

std::string
backtrace_string(int skip /*=0*/) {
	std::string bt_str;

	size_t const callstack_size = 128;
	void* callstack[callstack_size];
	int frames = backtrace(callstack, callstack_size);
	char** strs = backtrace_symbols(callstack, frames);
	for ( int i = skip; i < frames; ++i ) {
		bt_str += demangle(strs[i]);
		bt_str += '\n';
	}
	free(strs);
	return bt_str;
}



#else // _WIN32, etc.

#endif
