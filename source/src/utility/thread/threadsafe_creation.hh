// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_utility_thread_threadsafe_creation_HH
#define INCLUDED_utility_thread_threadsafe_creation_HH

// Unit headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Boost headers
#include <boost/function.hpp>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#include <thread>
#include <mutex>
#include <atomic>
#include <memory>
#endif


// C++ headers
#include <map>

namespace utility {
namespace thread {

// Macro for passing a mutex (which may or may not exist) to the
// safe creation functions
#if defined MULTI_THREADED
#define SAFELY_PASS_MUTEX(X) X
#ifdef OLDER_GCC
#define SAFELY_PASS_THREADSAFETY_BOOL(X) X
#else
#define SAFELY_PASS_THREADSAFETY_BOOL(X) false
#endif
#else
#define SAFELY_PASS_MUTEX(X) false
#define SAFELY_PASS_THREADSAFETY_BOOL(X) false
#endif

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
inline
void
safely_create_singleton(
	boost::function< T * () > creation_func,
#if defined MULTI_THREADED
	std::atomic< T * > & instance
#else
	T * & instance
#endif
) {

#ifdef MULTI_THREADED
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
	// not multithreaded; standard singleton instantiation logic
	if ( ! instance ) {
		instance = creation_func();
	}
#endif

}

/// @brief Safely instantiate a singleton class in a (possibly)
/// multithreaded context.  This version works with shared_ptrs.
///
/// @details In the non-multithreaded case, this simply checks the
/// singleton's instance member; in the multithreaded case, it
/// checks the instance member, then it obtains the singleton's
/// instance-creation mutex, then it checks the instance member
/// again, to ensure that no other thread has already created the
/// instance, it creates the instance, and then it releases the mutex.
///
/// @note This function contains ifdef'd special-case logic for GCC 4.8 and
/// 4.9, which lacked the C++11-standard functions std::atomic_load(std::shared_ptr) and
/// std::atomic_store(std::shared_ptr).  The workaround uses an std::atomic_bool for
/// data loading checks, and should result in very little cost in performance or memory.
/// (Thank you to Andrew Leaver-Fay for suggesting the workaround.)
///
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- modified from safely_create_singleton().
template < class T >
inline
void
safely_create_load_once_object_by_OP(
	boost::function< utility::pointer::shared_ptr< T > () > creation_func,
	utility::pointer::shared_ptr< T > & instance,
#if defined MULTI_THREADED
	std::mutex &mut,
#ifdef OLDER_GCC
	std::atomic_bool &atomicbool
#else
	bool const /*dummy2*/
#endif
#else
	bool const /*dummy*/,
	bool const /*dummy2*/
#endif
) {

#ifdef MULTI_THREADED
#ifdef OLDER_GCC
    std::atomic_bool local_atomicbool( atomicbool.load( std::memory_order_relaxed ) );
    std::atomic_thread_fence( std::memory_order_acquire );
    if( !local_atomicbool ) {
	std::lock_guard< std::mutex > lock(mut);
	local_atomicbool = atomicbool.load( std::memory_order_relaxed );
	if( !local_atomicbool ) {
	    utility::pointer::shared_ptr< T > local_instance( creation_func() );
	    instance = local_instance; //Should be okay, since protected by mutex -- no other read or write occurs without this mutex.
	    std::atomic_thread_fence( std::memory_order_release );
	    atomicbool.store( true, std::memory_order_relaxed );
	}
    }
#else // !OLDER_GCC
	// threadsafe version that uses c++11 interface
	// taken from http://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
	utility::pointer::shared_ptr< T > local_instance( std::atomic_load_explicit( &instance, std::memory_order_relaxed ) );
	std::atomic_thread_fence( std::memory_order_acquire );
	if ( ! local_instance ) {
		std::lock_guard< std::mutex > lock( mut );
		local_instance = std::atomic_load_explicit( &instance, std::memory_order_relaxed );
		if ( ! local_instance ) {
			local_instance = creation_func();
			std::atomic_thread_fence( std::memory_order_release );
			std::atomic_store_explicit( &instance, local_instance, std::memory_order_relaxed );
		}
	}
#endif // ifdef OLDER_GCC
#else
	// not multithreaded; standard singleton instantiation logic
	if ( ! instance ) {
		instance = creation_func();
	}
#endif

}


#ifdef MULTI_THREADED

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
template < class T, class L >
typename std::map< L, utility::pointer::shared_ptr< T > >::const_iterator
create_and_insert(
	typename boost::function< utility::pointer::shared_ptr< T > () > builder,
	utility::thread::ReadWriteMutex & wrm,
	L const & tname,
	typename std::map< L, utility::pointer::shared_ptr< T > > & tmap
)
{
	utility::pointer::shared_ptr< T > newT = builder();
	utility::thread::WriteLockGuard lock( wrm );
	typename std::map< L, utility::pointer::shared_ptr< T > >::const_iterator iter;

	iter = tmap.find( tname );
	if ( iter == tmap.end() ) {
		iter = tmap.insert( std::make_pair( tname, newT ) ).first;
	}
	return iter;
}

#endif

/// @brief Check for a string map key in a map of string->owning pointers.  If the key is not present, create an object by owning pointer
/// and insert it in the map with the given key.
/// @details Uses a utility::thread::ReadWriteMutex to allow many threads to read/check for the key,
/// and only one thread to write/create the object with the key.  Returns the object from the map.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
template< class T, class K >
typename utility::pointer::shared_ptr< T >
safely_check_map_for_key_and_insert_if_absent (
	typename boost::function< utility::pointer::shared_ptr< T > () > builder,
#ifdef MULTI_THREADED
	utility::thread::ReadWriteMutex & wrm,
#else
	bool const /*wrm*/,
#endif
	K const & tname,
	typename std::map< K, utility::pointer::shared_ptr< T > > & tmap
) {
	typename std::map< K, utility::pointer::shared_ptr< T > >::const_iterator iter, iter_end;
	{ //Scope for the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( wrm );
#endif
		iter = tmap.find(tname);
		iter_end = tmap.end();
	} //Release read lock
	if ( iter == iter_end ) {
#ifdef MULTI_THREADED
		iter = create_and_insert( builder, wrm, tname, tmap );
#else
		utility::pointer::shared_ptr< T > newobj( builder() );
		iter = tmap.insert( std::make_pair( tname, newobj ) ).first;
#endif
	}
	return iter->second;
}

} //namespace thread
} //namespace utility

#endif
