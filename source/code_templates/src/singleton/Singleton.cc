// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <--path--/--class--.hh>


// Unit headers

// Project header
#include <core/types.hh>

// Utility headers

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {


#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< --class-- >::singleton_mutex_ {};
template <> std::atomic< --class-- * > utility::SingletonBase< --class-- >::instance_( 0 );
#else
template <> --class-- * utility::SingletonBase< --class-- >::instance_( 0 );
#endif

}  // namespace utility


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );


--namespace--

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
--class--::--class--()
{

}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
--class-- *
--class--::create_singleton_instance()
{
	return new --class--;
}

--end_namespace--
