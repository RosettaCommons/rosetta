// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/thread/backwards_thread_local.hh
/// @brief  File to provide backwards compatibility for the thread_local keyword
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- modified for clang 3.7.0 case.


#ifdef MULTI_THREADED

// You need thread_local for multithreading, so if this causes a compiler error
// you'll need to update your compiler. No amount of preprocessor magic will help
#define THREAD_LOCAL thread_local

#else

// To avoid issues with compilers which don't support thread_local, don't use it
// if we don't need it.
#define THREAD_LOCAL

#endif
