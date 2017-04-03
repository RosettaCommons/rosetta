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
#include <utility/excn/Exceptions.hh>
#include <utility/CSI_Sequence.hh>

#include <iostream>

static bool throw_the_next_time_an_assertion_failure_is_hit( false );

void set_throw_on_next_assertion_failure()
{
	throw_the_next_time_an_assertion_failure_is_hit = true;
}

bool maybe_throw_on_next_assertion_failure( char const * condition )
{
	if ( throw_the_next_time_an_assertion_failure_is_hit ) {
		throw_the_next_time_an_assertion_failure_is_hit = false;
		throw utility::excn::EXCN_Msg_Exception( std::string( "assertion failure hit:" ) + condition );
	}
	return false;
}

bool
handle_assert_failure( char const * condition, std::string const & file, int const line ) {

	// Instead of printing a backtrace and halting the program, throw an exception.
	// Look, don't rely on this functionality in anything besides your unit tests.
	maybe_throw_on_next_assertion_failure( condition );

	std::cerr << utility::CSI_Reset() << utility::CSI_Red() << utility::CSI_Bold();
	std::cerr << "\nERROR: Assertion `" << condition << "` failed.\n";
	std::cerr << "ERROR:: Exit from: " << file << " line: " << line << "\n";
	std::cerr << utility::CSI_Reset();

	print_backtrace( condition );

	std::cerr << std::endl;

#ifdef __clang_analyzer__
	abort(); // To make the compiler happy on release-mode builds
#else
	return false;
#endif
}
