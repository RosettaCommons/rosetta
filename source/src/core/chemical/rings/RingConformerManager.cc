// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/rings/RingConformerManager.cc
/// @brief   Method definitions for RingConformerManager.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/rings/RingConformer.hh>
#include <core/chemical/rings/RingConformerManager.hh>
#include <core/chemical/rings/ring_conformer_io.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/database/open.hh>

// C++ headers
#include <map>
#include <sstream>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::chemical::rings::RingConformerManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< RingConformerManager >::singleton_mutex_ {};
template <> std::atomic< RingConformerManager * > utility::SingletonBase< RingConformerManager >::instance_( 0 );
#else
template <> RingConformerManager * utility::SingletonBase< RingConformerManager >::instance_( 0 );
#endif

}  // namespace utility


namespace core {
namespace chemical {
namespace rings {

using namespace core;


// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
utility::vector1< RingConformer > const &
RingConformerManager::conformers_for_ring_size( core::Size ring_size )
{
	return get_instance()->get_conformers_for_ring_size( ring_size );
}


// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RingConformerManager::RingConformerManager() {}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
RingConformerManager *
RingConformerManager::create_singleton_instance()
{
	return new RingConformerManager;
}


// Get the conformers requested, creating them if necessary.
// Called by the public static method conformers_for_ring_size().
utility::vector1< RingConformer > const &
RingConformerManager::get_conformers_for_ring_size( core::Size ring_size )
{
	using namespace std;
	using namespace utility;

	// Only create sets one time, as needed, for each ring size.
	if ( ! conformers_.count( ring_size ) ) {
		stringstream filename( stringstream::out );
		filename << "chemical/ring_conformer_sets/" << ring_size << "-membered_ring_conformers.data";
		vector1< RingConformer > conformers( read_conformers_from_database_file_for_ring_size(
			basic::database::full_name( filename.str() ), ring_size ) );
		conformers_.insert( make_pair( ring_size, conformers ) );
	}
	return conformers_[ ring_size ];
}

}  // namespace rings
}  // namespace chemical
}  // namespace core
