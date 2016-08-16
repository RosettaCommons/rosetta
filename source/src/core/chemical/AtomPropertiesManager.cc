// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/AtomPropertiesManager.cc
/// @brief   Method definitions for AtomPropertiesManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/chemical/AtomPropertiesManager.hh>
#include <core/chemical/AtomProperty.hh>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::chemical::AtomPropertiesManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< AtomPropertiesManager >::singleton_mutex_ {};
template <> std::atomic< AtomPropertiesManager * > utility::SingletonBase< AtomPropertiesManager >::instance_( 0 );
#else
template <> AtomPropertiesManager * utility::SingletonBase< AtomPropertiesManager >::instance_( 0 );
#endif

}  // namespace utility


namespace core {
namespace chemical {

// Public methods /////////////////////////////////////////////////////////////
AtomProperty const &
AtomPropertiesManager::property_from_string( std::string const & property )
{
	return get_instance()->string_to_property_map().find( property )->second;
}

std::string const &
AtomPropertiesManager::string_from_property( AtomProperty const property )
{
	return get_instance()->property_to_string_map().find( property )->second;
}

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
AtomPropertiesManager::AtomPropertiesManager() {}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
AtomPropertiesManager *
AtomPropertiesManager::create_singleton_instance()
{
	return new AtomPropertiesManager;
}

}  // namespace chemical
}  // namespace core
