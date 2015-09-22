// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/thread/backwards_thread_local.hh
/// @brief  File to provide backwards compatibility for the thread_local keyword
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- modified for clang 3.7.0 case.

#if (!defined CXX11 || !defined MULTI_THREADED) 
//Case 1: thread_local keyword not defined
#define THREAD_LOCAL

#else
//Case 2: thread_local keyword defined
#define THREAD_LOCAL thread_local

#endif
