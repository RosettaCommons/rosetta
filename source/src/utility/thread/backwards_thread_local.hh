// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/thread/backwards_thread_local.hh
/// @brief  File to provide backwards compatibility for the thread_local keyword
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#if !defined CXX11 || !defined MULTI_THREADED

// If we're not using both the cxx11 and multithreaded macros, then assume
// we need to erase the word "thread_local" from the code.  Maybe this
// is a bad idea?
#define thread_local

#endif
