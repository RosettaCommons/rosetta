// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/CSI_Sequence.cc
/// @brief  Terminal ASCII codes
/// @author Sergey Lyskov, @modified by Caleb Geniesse


#include <utility/CSI_Sequence.hh>

#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif

#include <cstdio>


namespace utility {

bool stdout_is_tty() {
#ifdef WIN32
	return _isatty(fileno(stdout));
#else
	return isatty(fileno(stdout));
#endif
}

// No mutexes, as this is really only going to change during setup
bool & CSI_Sequence::suppress_CSI_seq() {
	// This is a deliberately leaked pointer to a heap bool to insure that it lasts until the very end of the program,
	// as colored text can potentially be output as a result of object destructors during program tear-down.
	static auto * do_suppression( new bool(!stdout_is_tty()) ); // Initial value, may be reset
	return *do_suppression;
}

} // utility
