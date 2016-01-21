// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_utility_thread_threadsafe_creation_HH
#define INCLUDED_utility_thread_threadsafe_creation_HH

// Unit headers
#include <utility/thread/threadsafe_creation.fwd.hh>

// Boost headers
#include <boost/function.hpp>

#ifdef MULTI_THREADED
#ifdef CXX11

#include <utility/thread/ReadWriteMutex.hh>
#include <thread>

#endif
#endif

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <map>

namespace utility {
namespace thread {

/// @brief Safely instantiate a singleton class in a (possibly)
/// multithreaded context.
///
/// @details In the non-multithreaded case, this simply checks the
/// singleton's instance member; in the multithreaded case, it
/// checks the instance member, then it obtains the singleton's
/// instance-creation mutex, then it checks the instance member
/// again, to ensure that no other thread has already created the
/// instance, it creates the instance, and then it releases the mutex.
///
/// Requires that class T defines:
///   static std::mutex & singleton_mutex(),
template < class T >
void
safely_create_singleton(
	boost::function< T * () > creation_func,
#if defined MULTI_THREADED && defined CXX11
	std::atomic< T * > & instance
#else
	T * & instance
#endif
) {

#ifdef MULTI_THREADED
#ifdef CXX11
	// threadsafe version that uses c++11 interface
	// taken from http://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
	T * local_instance = instance.load( std::memory_order_relaxed );
	std::atomic_thread_fence( std::memory_order_acquire );
	if ( ! local_instance ) {
		std::lock_guard< std::mutex > lock( T::singleton_mutex() );
		local_instance = instance.load( std::memory_order_relaxed );
		if ( ! local_instance ) {
			local_instance = creation_func();
			std::atomic_thread_fence( std::memory_order_release );
			instance.store( local_instance, std::memory_order_relaxed );
		}
	}
#else
  // ok, multithreaded w/o cxx11
	// not actually threadsafe!
	// http://www.cs.umd.edu/~pugh/java/memoryModel/DoubleCheckedLocking.html
	if ( ! instance ) {
		instance = creation_func();
	}
#endif
#else
	// not multithreaded; standard singleton instantiation logic
	if ( ! instance ) {
		instance = creation_func();
	}
#endif

}


#ifdef MULTI_THREADED
#ifdef CXX11

/// @brief Safely create and insert an object into a map that is guarded by
/// a utility::thread::ReadWriteMutex.  Uses boost::function to invoke
/// the creation function only if needed (i.e. it will not invoke the creation
/// function if the item has already been created by the time the WriteLockGuard
/// has been acquired.
///
/// @details First create a WriteLockGuard for the input ReadWriteMutext, blocking
/// until all other threads stop reading from or writing to the tmap.  Then search
/// the tmap to make sure that object being created (and identified by "tname") has
/// not already been inserted into the map (after all, another thread may have called
/// this very function before this thread did, and may have already created the
/// desired object).  If the object is still not in the tmap, then it is safe to
/// create the object and to insert it into the map.
///
/// This function returns an iterator for the newly-inserted map element.
template < class T >
typename std::map< std::string, utility::pointer::shared_ptr< T > >::const_iterator
create_and_insert(
	typename boost::function< utility::pointer::shared_ptr< T > () > builder,
	utility::thread::ReadWriteMutex & wrm,
	std::string const & tname,
	typename std::map< std::string, utility::pointer::shared_ptr< T > > & tmap
)
{
	utility::thread::WriteLockGuard lock( wrm );
	typename std::map< std::string, utility::pointer::shared_ptr< T > >::const_iterator iter;

	iter = tmap.find( tname );
	if ( iter == tmap.end() ) {
		utility::pointer::shared_ptr< T > newT = builder();
		iter = tmap.insert( std::make_pair( tname, newT )).first;
	}
	return iter;
}

#endif
#endif


}
}

#endif
