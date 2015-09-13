// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfoManager.cc
/// @brief   Method definitions for CarbohydrateInfoManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit headers
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

// Utility header
#include <utility/exit.hh>

// Basic header
#include <basic/database/open.hh>

// C++ header
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::chemical::carbohydrates::CarbohydrateInfoManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< CarbohydrateInfoManager >::singleton_mutex_ {};
template <> std::atomic< CarbohydrateInfoManager * > utility::SingletonBase< CarbohydrateInfoManager >::instance_( 0 );
#else
template <> CarbohydrateInfoManager * utility::SingletonBase< CarbohydrateInfoManager >::instance_( 0 );
#endif

}  // namespace utility


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
bool
CarbohydrateInfoManager::is_valid_sugar_code( std::string const & code )
{
	return get_instance()->code_to_root_map().count( code );
}

std::string const &
CarbohydrateInfoManager::root_from_code( std::string const & code )
{
	if ( ! is_valid_sugar_code( code ) ) {
		utility_exit_with_message( "CarbohydrateInfoManager::root_from_code(): " +
			code + " is not a valid 3-letter code for carbohydrates.  " +
			"Has it been added to database/chemical/carbohydrates/codes_to_roots.map?");
	}
	return get_instance()->code_to_root_map().find( code )->second;
}

char
CarbohydrateInfoManager::ring_affix_from_ring_size( core::Size ring_size )
{
	return get_instance()->ring_size_to_morphemes_map().find( ring_size )->second.first;
}

std::string const &
CarbohydrateInfoManager::morpheme_from_ring_size( core::Size ring_size )
{
	return get_instance()->ring_size_to_morphemes_map().find( ring_size )->second.second;
}

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
CarbohydrateInfoManager::CarbohydrateInfoManager() {}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
CarbohydrateInfoManager *
CarbohydrateInfoManager::create_singleton_instance()
{
	return new CarbohydrateInfoManager;
}

// Get the map of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root requested,
// creating them if necessary.
// Called by the public static method root_from_code().
std::map< std::string, std::string > const &
CarbohydrateInfoManager::code_to_root_map()
{
	// Only creates map one time, as needed.
	if ( code_to_root_map_.empty() ) {
		code_to_root_map_ = read_codes_and_roots_from_database_file(
			basic::database::full_name( "chemical/carbohydrates/codes_to_roots.map" ) );
	}
	return code_to_root_map_;
}

// Get the map of carbohydrate ring sizes and their 1-letter affixes and morphemes requested, creating it if
// necessary.
// Called by the public static methods ring_affix_from_ring_size() and morpheme_from_ring_size().
std::map< core::Size, std::pair< char, std::string > > const &
CarbohydrateInfoManager::ring_size_to_morphemes_map()
{
	// Only creates map one time, as needed.
	if ( ring_size_to_morphemes_map_.empty() ) {
		ring_size_to_morphemes_map_ = read_ring_sizes_and_morphemes_fromt_database_file(
			basic::database::full_name( "chemical/carbohydrates/ring_size_to_morphemes.map" ) );
	}
	return ring_size_to_morphemes_map_;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
